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
import nzrec

import utils, params

pd.options.display.max_columns = 10

##########################################
### parcels

# way_id = 3081022
# way_id = 3135527
# f = booklet.open(params.sites_land_cover_reductions_path)


if __name__ == '__main__':
    ways = []
    with booklet.open(params.sites_land_cover_reductions_path) as f:
        for way_id in f:
            ways.append(way_id)

    # w0 = nzrec.Water(params.nzrec_data_path)

    # lake_segs = set(way_id for way_id, v in w0._way_tag.items() if v['Median flow'] is None)

    with concurrent.futures.ProcessPoolExecutor(max_workers=8, mp_context=mp.get_context("spawn")) as executor:
        with booklet.open(params.sites_land_cover_reductions_path) as f:
            # print(len(f))
            futures = {}
            for way_id in ways:
                data = f[way_id]
                f1 = executor.submit(utils.process_land_cover_per_catch, way_id, data)
                futures[f1] = way_id

        # Save results
        with booklet.open(params.sites_land_cover_catch_path, 'n', value_serializer='gpd_zstd', key_serializer='uint4', n_buckets=1607) as land_loads_dict:
            counter = 0
            for future in concurrent.futures.as_completed(futures):
                way_id = futures[future]
                run_result = future.result()
                counter += 1
                print(counter)
                land_loads_dict[way_id] = run_result

    # results0 = pd.concat(results_list)

    ## Save example catchment comparisons
    way_id = 3135527
    with booklet.open(params.sites_land_cover_reductions_path) as f:
        red = f[way_id]

    red.to_file(params.example_catch_land_cover_path)

    with booklet.open(params.sites_land_cover_catch_path) as f:
        red = f[way_id]

    red.to_file(params.example_catch_land_cover_rec_path)





##############################################
### Testing

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






























































































