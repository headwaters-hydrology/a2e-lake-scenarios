#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 08:00:37 2024

@author: mike
"""
import os
import geopandas as gpd
import pandas as pd
import numpy as np
import xarray as xr
import booklet

import utils, params

pd.options.display.max_columns = 10

####################################################
### Combine yields and reductions


def combine_red_yields():

    awm_red = pd.read_feather(params.awm_red_path).set_index('land_cover')
    awm_yields = pd.read_feather(params.awm_yields_path).set_index('land_cover')

    awm = pd.concat([awm_yields, awm_red], axis=1)
    awm['area_ha'] = 0

    combo_dict = {}
    combo_list = []
    with booklet.open(params.lakes_land_cover_reductions_path) as r:
        with booklet.open(params.lakes_yields_blt_path) as y:
            for lake_id, yields in y.items():

                if 'farm_type' in yields.columns:
                    yields = yields.drop('farm_type', axis=1)

                red = r[lake_id]

                # if lake_id == 53532:
                #     d

                combo0 = red.merge(yields.drop(['geometry', 'area_ha'], axis=1), on=['land_cover', 'typology'])
                combo0['land_use'] = combo0['land_cover'].replace(params.combine_land_covers)

                combo1 = combo0.copy()
                combo1['LFENZID'] = lake_id
                combo_list.append(combo1)

                combo2 = combo0[['land_use', 'area_ha', 'geometry']].dissolve('land_use', aggfunc='sum')
                values = utils.area_weighted_mean(combo0, 'land_use', ['phosphorus_yield', 'nitrogen_yield', 'phosphorus_reduction', 'nitrogen_reduction']).round(3)
                combo3 = pd.concat([combo2, values], axis=1)

                # Remove geometry for now, but I might need it later...
                combo4 = pd.DataFrame(combo3.drop('geometry', axis=1))
                combo5 = pd.concat([combo4, awm])
                combo5 = combo5.loc[~combo5.index.duplicated()]

                combo_dict[lake_id] = combo5.reindex(params.lu_order)

    ## Save results
    with booklet.open(params.lakes_lc_red_yields_path, 'n', value_serializer='pickle', key_serializer='uint4', n_buckets=1607) as f:
        for catch_id, data in combo_dict.items():
            f[catch_id] = data

    combo6 = pd.concat(combo_list).reset_index(drop=True)
    combo6.to_file(params.lakes_lc_red_yields_gpkg)


















































