#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 09:51:53 2024

@author: mike
"""
import os
from gistools import vector, rec
import geopandas as gpd
import pandas as pd
import numpy as np
from copy import copy
import nzrec
import booklet
import orjson
import geobuf
import shapely
from shapely.geometry import Point, Polygon, box, LineString, mapping

import utils, params

pd.options.display.max_columns = 10


#############################################
### Rivers


rec_tags = ['Climate class', 'Geology class', 'Landcover class', 'Network position class', 'Topography class', 'Valley landform class']


def lawa_wq_sites_for_app():
    """

    """
    wq_sites = gpd.read_file(params.wq_sites_raw_path)
    wq_data = pd.read_feather(params.wq_data_raw_path).replace({'Indicator': {'DIN': 'NO3N'}})
    sites_remove = pd.read_csv(params.sites_remove_path).LawaSiteID.values

    wq_sites = wq_sites[~wq_sites.LawaSiteID.isin(sites_remove)].copy()
    wq_data = wq_data[~wq_data.LawaSiteID.isin(sites_remove)].copy()

    ## Add in the source of flow from the REC
    w0 = nzrec.Water(params.nzrec_data_path)

    extra = []
    for seg in wq_sites.nzsegment:
        tags = w0._way_tag[seg]
        d1 = {key: val for key, val in tags.items() if key in rec_tags}
        d1['nzsegment'] = seg
        climate, topo, geo = tags['AmmendedCSOFG'].split('/')
        # if sof0.count('/') > 1:
        #     l = sof0.split('/')
        #     sof0 = l[0] + '/' + l[1]

        if tags['Climate class'] is not None:
            climate = tags['Climate class']

        if tags['Topography class'] is not None:
            topo = tags['Topography class']

        if tags['Geology class'] is not None:
            geo = tags['Geology class']

        d1['source'] = climate + '/' + topo

        if (climate in ('WD', 'CD')) and (geo in ('SS', 'VA', 'VB')):
            d1['Periphyton class'] = 'productive'
        else:
            d1['Periphyton class'] = 'default'

        extra.append(d1)

    extra_df = pd.DataFrame(extra)

    wq_sites = wq_sites.merge(extra_df, on='nzsegment').drop_duplicates('LawaSiteID')

    ## Only indicators we're interested in
    wq_data0 = wq_data[wq_data.Indicator.isin(params.indicators)].set_index(['LawaSiteID', 'Indicator'])

    ## Must have recent data
    grp = wq_data0.groupby(['LawaSiteID', 'Indicator'])

    max_dates = grp.SampleDateTime.max()
    recent_bool = (max_dates > '2021-07-01')
    # recent_bool.name = 'remove'

    max_dates = max_dates[recent_bool] + pd.DateOffset(seconds=1)
    max_dates.name = 'max_date'

    ## Must have at least 45 samples in the past 5 years
    min_dates = pd.to_datetime((max_dates - pd.DateOffset(years=5) - pd.offsets.MonthBegin()).dt.date)
    min_dates.name = 'min_date'
    min_max_dates = pd.concat([min_dates, max_dates], axis=1)

    data_list = []
    for i, data in min_max_dates.iterrows():
        data1 = wq_data0.loc[i]
        data2 = data1[(data1.SampleDateTime >= data.min_date) & (data1.SampleDateTime <= data.max_date)].reset_index()
        if len(data2) >= 45:
            data_list.append(data2)

    wq_data1 = pd.concat(data_list)

    ## Filter out sites
    sites = wq_data1['LawaSiteID'].unique()
    wq_sites1 = wq_sites[wq_sites.LawaSiteID.isin(sites)].copy()

    ## Make site name
    site_names_list = []
    for i, row in wq_sites1.iterrows():
        try:
            int_name = int(row['SiteID'])
            name = row['CouncilSiteID']
        except ValueError:
            name = row['SiteID']

        site_names_list.append(name)

    wq_sites1['site_name'] = site_names_list

    ## Correct censored data
    wq_data1.loc[wq_data1.Symbol == '<', 'Value (numeric)'] = wq_data1.loc[wq_data1.Symbol == '<', 'Value (numeric)'] * 0.5
    wq_data1.loc[wq_data1.Symbol == '>', 'Value (numeric)'] = wq_data1.loc[wq_data1.Symbol == '>', 'Value (numeric)'] * 1.5

    wq_data2 = wq_data1[['LawaSiteID', 'Indicator', 'SampleDateTime', 'Value (numeric)']].rename(columns={'Value (numeric)': 'value'}).groupby(['LawaSiteID', 'Indicator', 'SampleDateTime'])['value'].mean().reset_index()

    ## Save data
    # Sites
    wq_sites2 = wq_sites1.set_index('LawaSiteID', drop=False)
    wq_sites2['tooltip'] = wq_sites2['site_name']

    gbuf = geobuf.encode(wq_sites2.__geo_interface__)

    with open(params.wq_sites_gbuf_path, 'wb') as f:
        f.write(gbuf)

    wq_sites1.to_file(params.wq_sites_gpkg_path)

    # Sites/Indicators
    sites_ind = wq_data2[['LawaSiteID', 'Indicator']].drop_duplicates()
    with booklet.open(params.wq_sites_ind_path, 'n', key_serializer='str', value_serializer='pickle', n_buckets=1607) as f:
        for site_id, inds in sites_ind.groupby('LawaSiteID')['Indicator']:
            f[site_id] = inds.values.tolist()

    # WQ data
    with booklet.open(params.wq_data_path, 'n', key_serializer='str', value_serializer='pd_zstd', n_buckets=1607) as f:
        for site_id, data in wq_data2.groupby('LawaSiteID'):
            f[site_id] = data.reset_index(drop=True)

    # Sites names
    with booklet.open(params.wq_sites_names_path, 'n', key_serializer='str', value_serializer='str', n_buckets=1607) as f:
        for i, row in wq_sites1.iterrows():
            f[row['LawaSiteID']] = row['site_name']



def test_sample_size_filter():
    """

    """
    # wq_sites = gpd.read_file(params.wq_sites_raw_path)
    wq_data = pd.read_feather(params.wq_data_raw_path)

    wq_data_extra = wq_data[wq_data.Indicator.isin(['NO3N', 'DIN'])].copy()
    wq_data_extra['Indicator'] = 'DIN + NO3N'

    wq_data = pd.concat([wq_data, wq_data_extra])

    ## Only indicators we're interested in
    ind = ['DRP', 'ECOLI', 'NO3N', 'TN', 'TP', 'DIN + NO3N']
    wq_data0 = wq_data[wq_data.Indicator.isin(ind)].set_index(['LawaSiteID', 'Indicator'])

    ## Must have recent data
    grp = wq_data0.groupby(['LawaSiteID', 'Indicator'])

    max_dates = grp.SampleDateTime.max()
    recent_bool = (max_dates > '2021-07-01')
    # recent_bool.name = 'remove'

    max_dates = max_dates[recent_bool] + pd.DateOffset(seconds=1)
    max_dates.name = 'max_date'

    ## Must have at least 54 samples in the past 5 years
    min_dates = pd.to_datetime((max_dates - pd.DateOffset(years=5) - pd.offsets.MonthBegin()).dt.date)
    min_dates.name = 'min_date'
    min_max_dates = pd.concat([min_dates, max_dates], axis=1)

    data_list = []
    for i, data in min_max_dates.iterrows():
        data1 = wq_data0.loc[i]
        data2 = data1.loc[(data1.SampleDateTime >= data.min_date) & (data1.SampleDateTime <= data.max_date), 'SampleDateTime'].reset_index()
        data_list.append(data2)
        # if len(data2) > 54:
        #     data_list.append(data2)

    wq_data1 = pd.concat(data_list)
    # grp = wq_data1.groupby(['LawaSiteID', 'Indicator'])

    # results_list = []
    # for i, data in grp:
    #     d

    count1 = wq_data1.groupby(['LawaSiteID', 'Indicator'])['SampleDateTime'].count()

    results_list = []
    for thres in range(45, 61):
        # print(thres)
        count2 = count1[count1 >= thres]
        count3 = count2.groupby('Indicator').count()
        count3.name = thres
        results_list.append(count3.to_frame())

    results = pd.concat(results_list, axis=1).transpose()










































