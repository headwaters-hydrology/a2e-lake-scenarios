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

import utils, params

pd.options.display.max_columns = 10

#################################################
### Process loads

# way_id = 11018765


def process_inflow_loads_lakes():
    """

    """
    ## Yields per segment
    yields0 = pd.read_csv(params.rivers_loads_csv).set_index('nzsegment').rename(columns={'CurrentDRP': 'DRP', 'CurrentNNN': 'NNN', 'CurrentTN': 'TN',  'CurrentTP': 'TP', 'catAreaKM2': 'area_km2', 'Qmean': 'q_mean'})

    q_mean0 = yields0['q_mean'].copy()
    area_ha = (yields0['area_km2'] * 1000).copy()

    yields0 = yields0.drop(['q_mean', 'area_km2'], axis=1)

    yields_list = []
    for ind in yields0.columns:
        s = yields0[ind] * area_ha
        s.name = ind
        yields_list.append(s)

    yields1 = pd.concat(yields_list, axis=1)

    yield_segs = set(yields1.index)

    ## Combine yields with lake outputs
    lake_yields_dict = {}
    # missing_segs = set()
    with booklet.open(params.lake_outflows_blt) as f:
        for LFENZID, segs in f.items():
            # diff = segs.difference(yield_segs)
            # if diff:
            #     missing_segs.update(diff)
            both_segs = segs.intersection(yield_segs)
            y1 = yields1.loc[list(both_segs)].sum()
            a1 = area_ha.loc[list(both_segs)].sum()

            d1 = y1.round(2).to_dict()
            d1['area'] = round(a1)
            lake_yields_dict[LFENZID] = d1

    with booklet.open(params.lake_load_inflows, 'n', value_serializer='pickle', key_serializer='uint4', n_buckets=10007) as f:
        for LFENZID, res in lake_yields_dict.items():
            f[LFENZID] = res






























































