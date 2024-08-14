#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 09:46:22 2024

@author: mike
"""
import booklet
import pandas as pd
import geopandas as gpd
import utils, params


################################################
### Parameters


agency = 'Waikato Regional Council'
region = 'waikato'

waikato_path = params.data_path.joinpath('waikato')
waikato_cropland_csv = waikato_path.joinpath('waikato_cropland_areas_lcdb5.csv')
waikato_land_covers_csv = waikato_path.joinpath('waikato_land_cover_areas_lcdb5.csv')

###############################################
### % of cropland

sites = gpd.read_file(params.wq_sites_gpkg_path)
sites1 = sites[sites.Region == region].copy()

sites_segs = sites1.nzsegment.unique()

results_list = []
with booklet.open(params.sites_land_cover_reductions_path) as r:
    for seg in sites_segs:
        data = r[seg]
        area_sum1 = data.groupby('land_cover')['area_ha'].sum()
        tot_area = area_sum1.sum()
        area_ratio1 = ((area_sum1/tot_area) * 100).round(1)
        area_ratio1.name = 'area_perc'
        combo1 = pd.concat([area_sum1, area_ratio1], axis=1)
        combo1['nzsegment'] = seg
        results_list.append(combo1.reset_index())

results1 = pd.concat(results_list)

results2 = pd.merge(sites1[['nzsegment', 'site_name']], results1, on='nzsegment').drop_duplicates(['nzsegment', 'land_cover']).drop('nzsegment', axis=1)

results3 = results2[results2.land_cover.isin(['Orchard, Vineyard or Other Perennial Crop', 'Short-rotation Cropland'])].copy()

sum1 = results3.groupby('site_name')['area_perc'].sum()

sum_sites = sum1[sum1 >= 1].index

results4 = results2[results2.site_name.isin(sum_sites)]

results5 = results4[results4.site_name.isin(sites1.site_name)].sort_values(['site_name', 'land_cover'])

results6 = results2[results2.site_name.isin(sites1.site_name)].sort_values(['site_name', 'land_cover'])

results5.to_csv(waikato_cropland_csv, index=False)
results6.to_csv(waikato_land_covers_csv, index=False)




































































