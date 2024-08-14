#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:13:59 2024

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
### comparisons

way_id = 3081022
way_id = 3135527

reaches = booklet.open(params.sites_reach_mapping_path)

lc_loads_dict = booklet.open(params.sites_land_cover_catch_loads_path)
data0 = lc_loads_dict[way_id]

cols = [col for col in data0.columns if 'area' not in col]

res_list = []
for col in cols:
    data3 = data0[col]/(data0['area_m2']/10000)
    res_list.append(data3)

res = pd.concat(res_list, axis=1)
res.columns = cols

lc_loads = gpd.read_feather(params.sites_land_cover_catch_loads_feather_path).set_index('nzsegment')

site_reaches = np.asarray(list(reaches[way_id]))
lc_loads1 = lc_loads.loc[site_reaches]

lc_loads2 = lc_loads1.drop(['geometry'], axis=1).groupby(['land_cover']).sum()

cols = [col for col in lc_loads2.columns if 'area' not in col]

res_list = []
for col in cols:
    data3 = lc_loads2[col]/(lc_loads2['area_m2']/10000)
    res_list.append(data3)

res = pd.concat(res_list, axis=1)
res.columns = cols



ton0 = pd.read_csv(params.ton_tn_loads_csv_path).set_index('nzsegment')
ton0.columns.name = 'lc'
ton_reaches = ton0.index.values
ton_site_reaches = site_reaches[np.isin(site_reaches, ton_reaches)]

# ton1 = ton1.loc[way_id]

ton2 = ton0.loc[ton_site_reaches]

ton1 = ton0.rename(columns=params.ton_lc_map).groupby('lc', axis=1).sum()

ton2 = ton1.loc[ton_site_reaches]

ton3 = ton0.loc[way_id]














































