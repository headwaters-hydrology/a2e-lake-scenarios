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
from shapely import intersection, difference, intersects, symmetric_difference, make_valid
import booklet
import orjson
import geobuf
import concurrent.futures
import multiprocessing as mp

import utils, params

pd.options.display.max_columns = 10

##########################################
### parcels

# parcels = gpd.read_file(utils.parcels_path)
# parcels = parcels[['id', 'geometry']].copy()

# catch_id = 3081022
# catch = catches[catch_id]


# special_typo = {3076139: {
#     'recommended': {
#         'Warm/Low/Well/Moist': {'total phosphorus': 27, 'total nitrogen': 18},
#         'Cool/Low/Well/Moist': {'total phosphorus': 12, 'total nitrogen': 7},
#         'Cool/Low/Light/Moist': {'total phosphorus': 20, 'total nitrogen': 12},
#         'Warm/Low/Light/Moist': {'total phosphorus': 20, 'total nitrogen': 12},
#         'Warm/Moderate/Light/Moist': {'total phosphorus': 20, 'total nitrogen': 12},
#         'Warm/Moderate/Well/Moist': {'total phosphorus': 20, 'total nitrogen': 12},
#         'High Producing Exotic Grassland': {'total phosphorus': 20, 'total nitrogen': 12},
#                                   },
#     }
#     }



line_break = '<br />'
bold_start = '<b>'
bold_end = '</b>'

yield_cols = ['phosphorus_yield', 'nitrogen_yield']

# special_base_name = '{way_id}_{special_name}_land_reductions.gpkg'


def make_tooltip(x):
    """

    """
    tt = '<b>Typology</b><br />{typo}<br /><b>Land Cover</b><br />{lc}'.format(typo=x['typology'], lc=x['land_cover'])

    return tt


lcdb = gpd.read_feather(params.lcdb_red_path).drop(['nitrogen_yield',  'phosphorus_yield'], axis=1)
snb_dairy = gpd.read_feather(params.snb_dairy_red_path).drop(['nitrogen_yield',  'phosphorus_yield'], axis=1)
snb_dairy = snb_dairy.explode(index_parts=False).reset_index(drop=True)

lcdb_awm = utils.area_weighted_mean(lcdb.replace({'land_cover': params.combine_land_covers}), 'land_cover', ['nitrogen_reduction',  'phosphorus_reduction'])
snb_dairy_awm = utils.area_weighted_mean(snb_dairy, 'land_cover', ['nitrogen_reduction',  'phosphorus_reduction'])

awm = pd.concat([snb_dairy_awm, lcdb_awm])
awm = awm.loc[~awm.index.duplicated()]
awm.loc['Native Vegetation'] = awm.loc['Other']

awm.reset_index().to_feather(params.awm_red_path)


if __name__ == '__main__':
    with concurrent.futures.ProcessPoolExecutor(max_workers=10, mp_context=mp.get_context("spawn")) as executor:
        futures = {}
        with booklet.open(params.lakes_catch_major_path) as catches:
            for lake_id, catch in catches.items():
                # if lake_id == 41314:
                #     d
                poly = gpd.GeoSeries([catch], crs=4326).to_crs(2193).iloc[0]

                lcdb1 = lcdb.loc[lcdb.sindex.query(poly, predicate="intersects")].copy()
                snb_dairy1 = snb_dairy.loc[snb_dairy.sindex.query(poly, predicate="intersects")].copy()

                f1 = executor.submit(utils.calc_reductions, lake_id, poly, lcdb1, snb_dairy1)
                futures[f1] = lake_id

        # Save results
        with booklet.open(params.lakes_land_cover_reductions_path, 'n', value_serializer='gpd_zstd', key_serializer='uint4', n_buckets=1607) as f:
            counter = 0
            for future in concurrent.futures.as_completed(futures):
                lake_id = futures[future]
                run_result = future.result()
                counter += 1
                print(counter)
                if run_result is not None:
                    f[lake_id] = run_result


























































































