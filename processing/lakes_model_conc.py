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


def process_model_conc_lakes():
    """

    """
    ## source data
    conc0 = pd.read_csv(params.lake_model_conc_csv, usecols=['LFENZID', 'Attribute', 'value'])

    ## Extract the indicator and stat
    conc_list = []
    for i, row in conc0.iterrows():
        attr = row.Attribute
        attr1 = attr.split('_')
        ind = attr1[0]
        stat = attr1[-1]
        if stat == 'Median':
            if ind == 'NH4N':
                if len(attr1) == 2:
                    conc_list.append([row.LFENZID, ind, row.value])
            else:
                conc_list.append([row.LFENZID, ind, row.value])

    conc1 = pd.DataFrame(conc_list, columns=['LFENZID', 'indicator', 'value'])

    ## convert units and round
    conc1.loc[conc1.indicator.isin(['TN', 'TP', 'NH4N']), 'value'] = (conc1.loc[conc1.indicator.isin(['TN', 'TP', 'NH4N']), 'value'] * 1000).round(2)
    conc1.loc[~conc1.indicator.isin(['TN', 'TP', 'NH4N']), 'value'] = conc1.loc[~conc1.indicator.isin(['TN', 'TP', 'NH4N']), 'value'].round(2)

    ## Save data for app
    with booklet.open(params.lake_model_conc_blt, 'n', value_serializer='pickle', key_serializer='uint4', n_buckets=10007) as f:
        for LFENZID, res in conc1.groupby('LFENZID'):
            f[LFENZID] = res.set_index('indicator')['value'].to_dict()































































