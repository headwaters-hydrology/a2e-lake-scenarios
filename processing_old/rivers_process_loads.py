#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 17:23:53 2023

@author: mike
"""
import os
import xarray as xr
from gistools import vector, rec
import geopandas as gpd
import pandas as pd
import numpy as np
import hdf5tools
import booklet
import nzrec
from copy import copy

import utils, params

pd.options.display.max_columns = 10

#################################################
### Process loads


# segs = [3087421, 3087553, 3088075, 3088017, 3088076, 3095138]
# way_id = 11018765

tags = ['Climate class', 'Geology class', 'Topography class']


def process_loads_rec():
    """

    """
    ## yield set 1 - TN, TP, DRP, NNN
    conc0 = pd.read_csv(params.rivers_conc_csv_path1).set_index('nzsegment')
    cols = [col for col in conc0.columns if ('Median' not in col)]
    conc1 = conc0[cols].copy()
    cols = [col for col in conc1.columns if ('Current' in col)]
    for col in cols:
        conc1[col] = conc1[col]*conc1.catAreaKM2*100

    yields1 = conc1.drop('catAreaKM2', axis=1)
    yields1.columns = [col.split('Current')[1] if 'NNN' not in col else 'NO3N' for col in cols]

    ## Calc load diferences from upstream segments
    w0 = nzrec.Water(params.nzrec_data_path)

    way = {k: v for k, v in w0._way.items()}
    way_index = {k: v for k, v in w0._way_index.items()}
    node_way = {k: v for k, v in w0._node_way_index.items()}

    # reaches = booklet.open(params.sites_reach_mapping_path)

    old_segs = yields1.index
    all_segs = np.asarray(list(way.keys()))
    missing_segs = all_segs[~np.isin(all_segs, old_segs)]
    missing_df = pd.DataFrame(0, index=missing_segs, columns=yields1.columns)
    yields2 = pd.concat([yields1, missing_df])
    yields = yields2.to_dict('dict')

    ## Calc diff flows
    with booklet.open(params.sites_reach_mapping_path) as reaches:
        up_loads_dict = {}
        for param, flows in yields.items():
            up_flows_dict = {}
            for way_id in reaches:
                # if way_id not in extra_end_segs:
                down_ways = set([way_id])
                down_flow = flows[way_id]

                # new_ways = set(way_index_3rd_up[way_id])
                new_ways = utils.get_directly_upstream_ways(way_id, node_way, way, way_index)

                if (down_flow is None) or (down_flow == 0):
                    up_flows_dict[way_id] = 0
                else:
                    up_flows_list = [flows[f] for f in new_ways]
                    calc_bool = any([f is None for f in up_flows_list])
                    if calc_bool:
                        up_flows_dict[way_id] = 0
                    else:
                        up_flow = sum(up_flows_list)

                        diff_flow = round(down_flow - up_flow, 4)
                        if diff_flow < 0:
                            up_flows_dict[way_id] = diff_flow
                            # print(diff_flow)
                        else:
                            up_flows_dict[way_id] = diff_flow

                while new_ways:
                    down_ways.update(new_ways)
                    old_ways = copy(new_ways)
                    new_ways = set()
                    for new_way_id in old_ways:
                        up_ways = set(way_index[new_way_id]).difference(down_ways)
                        new_ways.update(up_ways)

                        down_flow = flows[new_way_id]
                        if (down_flow is None) or (down_flow == 0):
                            up_flows_dict[new_way_id] = 0
                        else:
                            up_flows_list = [flows[f] for f in up_ways]
                            calc_bool = any([f is None for f in up_flows_list])
                            if calc_bool:
                                up_flows_dict[new_way_id] = 0
                            else:
                                up_flow = sum(up_flows_list)

                                diff_flow = round(down_flow - up_flow, 4)
                                if diff_flow <= 0:
                                    up_flows_dict[new_way_id] = diff_flow
                                    # print(diff_flow)
                                else:
                                    up_flows_dict[new_way_id] = diff_flow

            up_loads_dict[param] = up_flows_dict

    up_loads_dict2 = {}
    for param, vals in up_loads_dict.items():
        for way_id, val in vals.items():
            if way_id in up_loads_dict2:
                up_loads_dict2[way_id].update({param: val})
            else:
                up_loads_dict2[way_id] = {param: val}

    with booklet.open(params.river_loads_rec_diff_path, 'n', key_serializer='uint4', value_serializer='orjson') as f:
        for way_id, vals in up_loads_dict2.items():
            f[way_id] = vals

    ## Save as csv
    up_loads_df = pd.DataFrame.from_dict(up_loads_dict2, 'index')
    up_loads_df.index.name = 'nzsegment'
    up_loads_df.to_csv(params.river_loads_rec_diff_csv_path)

    parameters = list(up_loads_df.columns)

    ## Calc the total load per site catch
    loads = {n: 0 for n in parameters}

    with booklet.open(params.site_catch_rec_loads_path, 'n', value_serializer='orjson', key_serializer='uint4', n_buckets=1001) as f:
        with booklet.open(params.sites_reach_mapping_path) as reaches:
            for catch_id, reaches in reaches.items():
                load = loads.copy()
                for reach in reaches:
                    vals = up_loads_dict2[reach]
                    for n, val in vals.items():
                        load[n] += val
                f[catch_id] = load







### Testing
# combo2 = pd.merge(combo1, ref_load1, on='nzsegment')

# res_dict = {}
# for col in combo1.columns:
#     c1 = combo2[col+'_x']/combo2[col+'_y']
#     c2 = (c1 < 1).sum()
#     res_dict[col] = c2

# catch_loads = booklet.open(params.site_catch_rec_loads_path)

# for catch_id, loads0 in catch_loads.items():
#     for param, load0 in loads0.items():
#         if load0 < 0:
#             raise ValueError(f'{catch_id}')

# catch_loads.close()





















































