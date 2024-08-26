#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 18:15:09 2024

@author: mike
"""
import os
import sys
import pathlib
import geopandas as gpd
import pandas as pd
import numpy as np
import booklet

pd.options.display.max_columns = 10



#########################################################
### Parameters

script_path = pathlib.Path(os.path.realpath(os.path.dirname(__file__)))

processing_path = script_path.parent.joinpath('processing')

sys.path.append(str(processing_path))

import params, utils

######################################################
### Combine for Ton

yields_list = []
with booklet.open(params.lakes_yields_blt_path) as f:
    for catch_id, data in f.items():
        data['LFENZID'] = catch_id
        yields_list.append(data)

yields0 = pd.concat(yields_list)

if 'farm_type' in yields0.columns:
    yields0 = yields0.drop('farm_type', axis=1)

red_list = []
with booklet.open(params.lakes_land_cover_reductions_path) as f:
    for catch_id, data in f.items():
        data['LFENZID'] = catch_id
        red_list.append(data)

red0 = pd.concat(red_list)

yields1 = yields0.drop('geometry', axis=1)
red1 = red0.drop(['geometry', 'area_ha'], axis=1)
combo1 = pd.merge(yields1, red1, on=['LFENZID', 'typology', 'land_cover'], how='inner')
combo1['new_category'] = combo1['land_cover'].replace(params.combine_land_covers)

# sites_data = gpd.read_file(params.wq_sites_gpkg_path)

# combo2 = pd.merge(sites_data.drop('geometry', axis=1), combo1, on='nzsegment').drop_duplicates(['nzsegment', 'typology'])
combo1.set_index(['LFENZID', 'typology', 'land_cover', 'new_category']).to_csv(params.lake_lc_data_csv)

## Land cover summary
combo3 = combo1.drop(['area_ha', 'LFENZID'], axis=1)

grp1 = combo3.groupby('typology')

values = grp1.mean(numeric_only=True).round(3)
lcs = grp1['land_cover'].first()

combo4 = pd.concat([lcs, values], axis=1).reset_index()
combo4['new_category'] = combo4['land_cover'].replace(params.combine_land_covers)

combo5 = combo4.set_index(['typology', 'land_cover', 'new_category'])

combo5.to_csv(params.lake_lc_summ_csv)






































