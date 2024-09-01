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
from shapely import intersection
import xarray as xr
import booklet
import pickle
import concurrent.futures
import multiprocessing as mp

import utils, params

pd.options.display.max_columns = 10

##########################################
### land use/cover


# features_cols = ['Climate', 'Slope', 'Drainage', 'Wetness']


# yields0 = pd.read_csv(params.yields_csv_path).rename(columns={'TN': 'nitrogen_yield', 'TP': 'phosphorus_yield'})
# slope_geo0 = gpd.read_file(params.yields_gpkg_path, layer='Slope', include_fields=['Srinivasan.slope.type']).rename(columns={'Srinivasan.slope.type': 'slope'})
# slope_geo0['geometry'] = slope_geo0['geometry'].simplify(10)

# moisture_geo0 = gpd.read_file(params.yields_gpkg_path, layer='Moisture', include_fields=['Srinivasan.Moisture.Type']).rename(columns={'Srinivasan.Moisture.Type': 'moisture'})
# moisture_geo0['geometry'] = moisture_geo0['geometry'].simplify(10)

### SnB
# snb_dairy0 = gpd.read_feather(params.snb_dairy_red_path).drop(['phosphorus_yield', 'nitrogen_yield', 'phosphorus_reduction', 'nitrogen_reduction'], axis=1)

# lcdb0 = gpd.read_feather(params.lcdb_clean_path)
# lcdb_yields = pd.read_csv(params.lcdb_yields_csv_path).drop(['sediment_yield'], axis=1)
# lcdb_yields = lcdb_yields[lcdb_yields.typology != 'Native Vegetation'].copy()
# lcdb2 = lcdb0.merge(lcdb_yields, on='typology', how='outer').reset_index(drop=True)
# lcdb2.loc[lcdb2['nitrogen_yield'].isnull(), ['nitrogen_yield', 'phosphorus_yield']] = 0

# lcdb_awm = utils.area_weighted_mean(lcdb2.replace({'land_cover': params.combine_land_covers}), 'land_cover', ['nitrogen_yield',  'phosphorus_yield'])

# # awm = pd.concat([snb_dairy_awm, lcdb_awm])
# # awm = awm.loc[~awm.index.duplicated()]
# awm = lcdb_awm.copy()
# awm.loc['Native Vegetation'] = awm.loc['Other']
# awm.loc['Native Vegetation', ['nitrogen_yield', 'phosphorus_yield']] = [2, 0.3]

# awm.reset_index().to_feather(params.awm_yields_path)


if __name__ == '__main__':
    with concurrent.futures.ProcessPoolExecutor(max_workers=10, mp_context=mp.get_context("spawn")) as executor:
        futures = {}
        with booklet.open(params.lakes_catch_major_path) as catches:
            for lake_id, catch in catches.items():
                # if lake_id == 53532:
                #     d
                poly = gpd.GeoSeries([catch], crs=4326).to_crs(2193).iloc[0]

                f1 = executor.submit(utils.combine_yields, lake_id, poly)
                futures[f1] = lake_id

        # Save results
        with booklet.open(params.lakes_yields_blt_path, 'n', value_serializer='gpd_zstd', key_serializer='uint4', n_buckets=1607) as f:
            counter = 0
            for future in concurrent.futures.as_completed(futures):
                lake_id = futures[future]
                error = future.exception()
                if error is not None:
                    print(lake_id)
                    executor.shutdown(wait=True, cancel_futures=True)
                    raise error

                run_result = future.result()
                counter += 1
                print(counter)
                if run_result is not None:
                    f[lake_id] = run_result

































































































