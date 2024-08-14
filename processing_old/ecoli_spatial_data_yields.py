#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 13:11:07 2022

@author: mike
"""
import os
from gistools import vector, rec
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely import intersection, Point
import xarray as xr
import booklet
import concurrent.futures
import multiprocessing as mp

import utils, params

pd.options.display.max_columns = 10

##########################################
### land use/cover

drain_types = ('PoorlyDrained', 'WellDrained')
elev_types = ('LowElev', 'HighElev')


def ecoli_yields():
    #### Prep data
    
    # lu_map = pd.read_csv(params.ecoli_lu_csv_path)
    
    calcs = pd.read_csv(params.ecoli_input_params).drop(['Q0.05', 'Q0.95'], axis=1)
    land_type_list = []
    for i, row in calcs.iterrows():
        split = row['land_type'].split('_')
        lu = split[0]
        len1 = len(split)
        if len1 == 2:
            if lu == 'Dairy':
                drain = split[1]
                for elev in elev_types:
                    land_type_list.append({'lu': lu, 'drain': drain, 'elev': elev, 'ecoli_factor': row['median']})
            elif lu == 'Native':
                elev = split[1]
                for drain in drain_types:
                    land_type_list.append({'lu': lu, 'drain': drain, 'elev': elev, 'ecoli_factor': row['median']})
    
        elif len1 == 3:
            elev = split[1]
            drain = split[2]
            land_type_list.append({'lu': lu, 'drain': drain, 'elev': elev, 'ecoli_factor': row['median']})
        else:
            for elev in elev_types:
                for drain in drain_types:
                    land_type_list.append({'lu': lu, 'drain': drain, 'elev': elev, 'ecoli_factor': row['median']})
    
    ecoli_factors0 = pd.DataFrame(land_type_list)

    # crop1 = ecoli_factors0[ecoli_factors0.lu == 'Cropland'].copy()
    # crop1['lu'] = 'Perennial Crop'
    # crop2 = crop1.copy()
    # crop2['lu'] = 'Short-rotation Crop'

    # ecoli_factors = pd.concat([ecoli_factors0[ecoli_factors0.lu != 'Cropland'], crop1, crop2])

    ecoli_factors = ecoli_factors0.replace({'lu': {'Cropland': 'Short-rotation Crop', 'OrchardVineyard': 'Perennial Crop', 'Native': 'Native Vegetation', 'SheepBeef': 'Sheep and Beef'}})

    ### Prep
    ## Elevation
    elev0 = xr.load_dataset(params.ecoli_elev_tif_path, engine='rasterio').squeeze()['band_data'].round().drop(['band', 'spatial_ref'])
    elev0['x'] = elev0['x'].astype('int32')
    elev0['y'] = elev0['y'].astype('int32')
    elev0.name = 'elev'
    
    ## Drainage
    drain0 = xr.load_dataset(params.ecoli_drain_tif_path, engine='rasterio').squeeze()['band_data'].drop(['band', 'spatial_ref'])
    drain0['x'] = drain0['x'].astype('int32')
    drain0['y'] = drain0['y'].astype('int32')
    drain0.name = 'drain'

    ## land use
    # lu0 = xr.load_dataset(params.ecoli_lu_tif_path, engine='rasterio').squeeze()['band_data'].round().drop(['band', 'spatial_ref'])
    # lu0['x'] = lu0['x'].astype('int32')
    # lu0['y'] = lu0['y'].astype('int32')
    # lu0.name = 'lu_id'
    # lu0 = xr.where(lu0 == 10, np.nan, lu0)
    
    ## Combo
    # combo0 = xr.merge([elev0, drain0, lu0])
    combo0 = xr.merge([elev0, drain0])
    
    combo1 = combo0.to_dataframe().dropna(how='all')
    combo1.loc[combo1.elev.isnull(), 'elev'] = 3
    # combo1.loc[combo1.lu_id.isnull(), 'lu_id'] = 3
    combo1.loc[combo1.drain.isnull(), 'drain'] = 3
    
    combo1['elev'] = combo1['elev'].astype('int32')
    combo1['drain'] = combo1['drain'].astype('int8')
    # combo1['lu_id'] = combo1['lu_id'].astype('int8')
    
    # combo2 = vector.xy_to_gpd(['elev', 'lu_id', 'drain'], 'x', 'y', combo1.reset_index(), 2193).to_crs(2193)
    combo2 = vector.xy_to_gpd(['elev', 'drain'], 'x', 'y', combo1.reset_index(), 2193).to_crs(2193)
    combo2.sindex

    resutls_dict = {}
    counter = 0
    with booklet.open(params.sites_yields_blt_path) as y:
        print(f'Number of sites: {len(y)}')
        for catch_id, data in y.items():
            counter += 1
            print(counter)
            results = utils.ecoli_spatial_yields(data, combo2, ecoli_factors)
            resutls_dict[catch_id] = results

    with booklet.open(params.ecoli_catch_yields_path, 'n', value_serializer='pickle', key_serializer='uint4', n_buckets=1607) as f:
        for catch_id, data in resutls_dict.items():
            f[catch_id] = data




# if __name__ == '__main__':
#     with concurrent.futures.ProcessPoolExecutor(max_workers=10, mp_context=mp.get_context("spawn")) as executor:
#         futures = {}
#         with booklet.open(params.sites_yields_blt_path) as y:
#             print(f'Number of sites: {len(y)}')
#             for catch_id in y:
#                 f1 = executor.submit(utils.ecoli_spatial_yields, catch_id, combo2, lu_map, ecoli_factors)
#                 futures[f1] = catch_id

#         # Save results
#         with booklet.open(params.ecoli_catch_yields_path, 'n', value_serializer='pd_zstd', key_serializer='uint4', n_buckets=1607) as f:
#             counter = 0
#             for future in concurrent.futures.as_completed(futures):
#                 way_id = futures[future]
#                 run_result = future.result()
#                 counter += 1
#                 print(counter)
#                 if run_result is not None:
#                     f[way_id] = run_result




def testing(results):
    """

    """
    results1 = results.reset_index().copy()

    # means1 = utils.area_weighted_mean(results, 'land_use', 'ecoli_yield')
    areas1 = results1.groupby('land_use').area_ha.sum()
    areas_tot1 = areas1.sum()
    areas_prop_init1 = areas1/areas_tot1
    areas_prop_scen1 = areas_prop_init1.copy()
    areas_prop_scen1['Dairy'] -= 0.1
    areas_prop_scen1['Sheep and Beef'] -= 0.15
    areas_prop_scen1['Exotic Forest'] += 0.20
    areas_prop_scen1['Native Vegetation'] += 0.05
    areas_diff1 = (areas_prop_scen1 - areas_prop_init1).round(2)
    areas_diff1.name = 'area_diff_ratio'
    areas_lost1 = areas_diff1[areas_diff1 < 0]
    areas_lost2 = results1[results1.land_use.isin(areas_lost1.index)]
    areas_lost2 = pd.merge(results1, areas_lost1, on='land_use')
    areas_lost2['area_diff'] = areas_lost2.area_diff_ratio * areas_lost2.area_ha * -1
    areas_lost2['area_ha'] = areas_lost2['area_ha'] - areas_lost2['area_diff']
    extra_areas1 = areas_lost2.groupby(['elev', 'drain'])['area_diff'].sum().reset_index()
    areas_gained1 = areas_diff1[areas_diff1 > 0]
    tot_gained = areas_gained1.sum()
    areas_gained_ratio = areas_gained1/tot_gained

    gained_list = []
    for lu, ratio in areas_gained_ratio.items():
        agr = extra_areas1.copy()
        agr['land_use'] = lu
        agr['area_diff_ratio'] = ratio
        gained_list.append(agr)

    areas_gained2 = pd.concat(gained_list)
    areas_gained2['extra_area'] = areas_gained2['area_diff'] * areas_gained2['area_diff_ratio']

    areas_gained3 = results1[results1.land_use.isin(areas_gained1.index)]
    areas_gained4 = pd.merge(areas_gained3, areas_gained2[['land_use', 'elev', 'drain', 'extra_area']], on=['land_use', 'elev', 'drain'])
    areas_gained4['area_ha'] = areas_gained4['area_ha'] + areas_gained4['extra_area']

    areas_no_change1 = areas_diff1[areas_diff1 == 0]
    areas_no_change2 = results1[results1.land_use.isin(areas_no_change1.index)]

    ## Combine
    combo1 = pd.concat([areas_no_change2, areas_lost2[['land_use', 'elev', 'drain', 'area_ha', 'ecoli_factor']], areas_gained4.drop('extra_area', axis=1)])




































































































