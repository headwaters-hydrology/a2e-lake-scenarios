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

indicators = ('TN', 'TP', 'CHLA', 'Secchi')


def combine_model_monitoring_data():

    combo_dict = {}
    with booklet.open(params.lake_model_conc_blt) as model:
        with booklet.open(params.lake_moni_conc_blt) as moni:
            for lake_id, data in model.items():
                data['measured'] = False
                if lake_id in moni:
                    moni_data = moni[lake_id]
                    moni_bool = True
                    for ind in indicators:
                        if ind not in moni_data:
                            moni_bool = False
                    if moni_bool:
                        data = moni_data
                        data['measured'] = True

                combo_dict[lake_id] = data

    ## Save results
    with booklet.open(params.lake_model_moni_conc_blt, 'n', value_serializer='orjson', key_serializer='uint2', n_buckets=1607) as f:
        for catch_id, data in combo_dict.items():
            f[catch_id] = data



















































