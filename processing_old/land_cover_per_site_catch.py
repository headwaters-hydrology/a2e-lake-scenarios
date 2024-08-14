#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 13:11:07 2022

@author: mike
"""
import os
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely import intersection, difference, intersects, symmetric_difference
import booklet
import orjson
import geobuf
import concurrent.futures
import multiprocessing as mp

import utils, params

pd.options.display.max_columns = 10

##########################################
### parcels

# way_id = 3135527
# catch = catches[way_id]


def land_cover_per_site_catch():
    # lc_loads_dict = {}
    lc_loads_list = []
    with booklet.open(params.sites_land_cover_catch_loads_path, 'n', value_serializer='pickle', key_serializer='uint4', n_buckets=1607) as lc_loads_dict:
        with booklet.open(params.river_loads_rec_diff_path) as loads:
            with booklet.open(params.sites_land_cover_catch_path) as f:
                for catch_id, data in f.items():
                    segs = data.nzsegment.values
                    seg_loads = {seg: loads[seg] for seg in segs}
                    seg_loads_df = pd.DataFrame.from_dict(seg_loads, 'index')
                    seg_loads_df.index.name = 'nzsegment'
                    data1 = data[['nzsegment', 'land_cover', 'geometry']].merge(seg_loads_df, on='nzsegment')
                    data1['area_m2'] = data1.geometry.area.round().astype(int)
                    lc_loads_list.append(data1)

                    data2 = data1.drop(['nzsegment', 'geometry'], axis=1).groupby('land_cover').sum()
                    data2[data2 < 0] = 0
                    lc_loads_dict[catch_id] = data2

    lc_loads = pd.concat(lc_loads_list).drop_duplicates(subset=['nzsegment']).copy()
    utils.gpd_to_feather(lc_loads, params.sites_land_cover_catch_loads_feather_path)
    lc_loads.to_file(params.sites_land_cover_catch_loads_gpkg_path)

    # results0 = pd.concat(results_list)

    ## Save example catchment comparisons
    # way_id = 3135527
    # with booklet.open(params.sites_land_cover_reductions_path) as f:
    #     red = f[way_id]

    # red.to_file(params.example_catch_land_cover_path)

    # with booklet.open(params.sites_land_cover_catch_path) as f:
    #     red = f[way_id]

    # red.to_file(params.example_catch_land_cover_rec_path)





##############################################
### Testing


def testing():
    # land_cover_dict = shelflet.open(utils.catch_lc_path, 'r')
    
    
    # with shelve.open('/media/nvme1/data/OLW/web_app/output/shelve_test.shelf') as t:
    #     for seg in land_cover_dict:
    #         lc2 = land_cover_dict[seg]
    #         t[seg] = lc2
    #         t.sync()
    
    
    
    # db = booklet.open(utils.catch_lc_path)
    
    # keys = list(db.keys())
    
    # db[9259625]
    
    
    # with booklet.open(utils.catch_lc_path) as lc_dict:
    #     with booklet.open(utils.catch_lc_pbf_path, 'n', value_serializer=None, key_serializer='uint4', n_buckets=1600) as lc_gbuf:
    #         for i, lc2 in lc_dict.items():
    #             gdf = lc2.to_crs(4326)
    #             gjson = orjson.loads(gdf.to_json())
    #             gbuf = geobuf.encode(gjson)
    #             lc_gbuf[i] = gbuf
    
    way_id = 3135527
    
    lc_loads_dict = booklet.open(params.sites_land_cover_catch_loads_path)
    data0 = lc_loads_dict[way_id]
    
    cols = [col for col in data0.columns if 'area' not in col]
    
    res_list = []
    for col in cols:
        data3 = data0[col]/(data0['area_m2']/10000)
        res_list.append(data3)
    
    res = pd.concat(res_list, axis=1)
    res.columns = cols
    
    res_list = list(land_cover_loads.values())
    res = pd.concat(res_list, axis=0)
    
    res2 = res.groupby('land_cover').sum()
    
    cols = [col for col in res2.columns if 'area' not in col]
    
    res_list = []
    for col in cols:
        data3 = res2[col]/(res2['area_m2']/10000)
        res_list.append(data3)
    
    res = pd.concat(res_list, axis=1)
    res.columns = cols
    res['area_km2'] = res2['area_m2']*0.000001
    
    res = res.round(2)
    
    res.to_csv(params.land_cover_loads_agg_csv_path)
    

    lc_loads = gpd.read_feather(params.sites_land_cover_catch_loads_feather_path)
    
    lc_loads1 = lc_loads.drop(['nzsegment', 'geometry'], axis=1).groupby(['land_cover']).sum()
    
    cols = [col for col in lc_loads1.columns if 'area' not in col]
    
    res_list = []
    for col in cols:
        data3 = lc_loads1[col]/(lc_loads1['area_m2']/10000)
        res_list.append(data3)
    
    res = pd.concat(res_list, axis=1)
    res.columns = cols
    res['area_km2'] = lc_loads1['area_m2']*0.000001
    
    res = res.round(2)
    
    res.to_csv(params.land_cover_loads_agg_csv_path)
    
    
    
    ton0 = pd.read_csv(params.ton_tn_loads_csv_path).set_index('nzsegment')
    
    ton1 = ton0.loc[way_id]





























































