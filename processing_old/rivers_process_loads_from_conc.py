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
    ## Determine all the segs
    # catches_minor = booklet.open(params.sites_catch_minor_path)
    # segs = np.asarray(list(catches_minor.keys()))
    # catches_minor.close()

    ## Import flows for load calcs
    flows_list = []
    append = flows_list.append
    with booklet.open(params.river_flows_rec_diff_path) as f:
        for catch_id, flow in f.items():
            append((catch_id, flow))

    flows = pd.DataFrame(flows_list, columns=['nzsegment', 'flow'])

    ## Conc set 1 - TN, TP, DRP, NNN
    conc0 = pd.read_csv(params.rivers_conc_csv_path1).set_index('nzsegment')
    cols = [col for col in conc0.columns if ('Median' in col)]
    conc1 = conc0[cols].copy()
    cols1 = [col[7:].split('Median')[0] for col in cols]
    index = cols1.index('NNN')
    cols1.remove('NNN')
    cols1.insert(index, 'NO3N')
    conc1.columns = cols1

    loads1 = pd.merge(flows, conc1, on='nzsegment', how='left')
    loads1.loc[loads1[cols1[0]].isnull(), cols1] = 0
    for col in cols1:
        loads1[col] = loads1[col] * loads1['flow']

    loads1 = loads1.drop('flow', axis=1).set_index('nzsegment')

    ## Calc set 2 - e.coli
    conc0 = pd.read_csv(params.rivers_conc_csv_path2).set_index('nzsegment')
    cols1 = ['ECOLI']
    conc1 = conc0[['CurrentQ50']].copy()
    conc1.columns = cols1
    loads2 = pd.merge(flows, conc1, on='nzsegment', how='left')
    loads2.loc[loads2[cols1[0]].isnull(), cols1] = 0
    for col in cols1:
        loads2[col] = loads2[col] * loads2['flow']

    loads2 = loads2.drop('flow', axis=1).set_index('nzsegment')

    ## Calc set 3 - clarity and turbidity
    conc0 = pd.read_csv(params.rivers_conc_csv_path3, usecols=['nzsegment', 'cumArea', 'CurrCor_cu'])
    conc0['SS'] = conc0.CurrCor_cu/conc0.cumArea
    loads3 = pd.merge(flows, conc0[['nzsegment', 'SS']], on='nzsegment', how='left').set_index('nzsegment')
    loads3.loc[loads3.SS.isnull(), 'SS'] = 0
    # loads3['sediment'] = loads3['sediment'] * loads3['flow']

    loads3 = loads3.drop('flow', axis=1)

    ## Combine
    combo1 = pd.concat([loads1, loads2, loads3], axis=1).round(4)

    ## Convert to web app parameter
    # cols2 = combo1.columns
    # for param, col in utils.indicator_dict.items():
    #     combo1[param] = combo1[col].copy()

    # combo1 = combo1.drop(cols2, axis=1)

    # combo1.columns = [col+'_current_load' for col in combo1.columns]

    ## Save as csv
    combo1.to_csv(params.river_loads_rec_diff_csv_path, compression='infer')

    ## Save as dict
    loads_dict = combo1.to_dict('index')

    with booklet.open(params.river_loads_rec_diff_path, 'n', value_serializer='orjson', key_serializer='uint4', n_buckets=500001) as f:
        for seg, val in loads_dict.items():
            f[seg] = val

    ## Calc the total load per site catch
    loads = {n: 0 for n in val}

    with booklet.open(params.site_catch_rec_loads_path, 'n', value_serializer='orjson', key_serializer='uint4', n_buckets=1001) as f:
        with booklet.open(params.sites_reach_mapping_path) as sites_reaches:
            for catch_id, reaches in sites_reaches.items():
                load = loads.copy()
                for reach in reaches:
                    vals = loads_dict[reach]
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



























































