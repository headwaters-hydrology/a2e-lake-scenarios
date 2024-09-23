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

indicators = ('TN', 'TP', 'CHLA', 'Secchi')

def process_achievable_concs():
    """

    """
    ## Calc N and P inflow ratios
    yields_ratios = {}
    with booklet.open(params.lakes_lc_red_yields_path) as f:
        for lake_id, data in f.items():
            # if lake_id == 47833:
            #     d
            area_ha = data['area_ha']['area_ha']
            tot_area = area_ha.sum()
            native_area = tot_area - area_ha['Other']
            yields_current = data['yield'].multiply(area_ha, axis=0).sum()
            yields_scenario = data.loc['Native Vegetation', 'yield'].multiply(native_area, axis=0)

            yields_ratio = yields_scenario/yields_current
            best_ratio = yields_ratio.round(6).to_dict()

            yields_scenario = data.loc['Dairy', 'yield'].multiply(native_area, axis=0)
            yields_ratio = yields_scenario/yields_current
            worst_ratio = yields_ratio.round(6).to_dict()

            yields_ratios[lake_id] = {'best': best_ratio, 'worst': worst_ratio}

    ## Estimate best/worst achievable scenario
    achieve_conc = {}
    with booklet.open(params.lake_model_moni_conc_blt) as f:
        with booklet.open(params.lake_metadata_blt) as meta:
            for lake_id, lake_meta in meta.items():
                if lake_id in f:
                    data = f[lake_id]

                    region = lake_meta['regional_council']
                    t = lake_meta['residence_time']
                    z_max = lake_meta['max_depth']

                    yields_ratio = yields_ratios[lake_id]
                    best_ratio = yields_ratio['best']
                    worst_ratio = yields_ratio['worst']

                    c_best = utils.est_ind_scenario_conc(data, best_ratio, region, t, z_max)
                    c_worst = utils.est_ind_scenario_conc(data, worst_ratio, region, t, z_max)
                    achieve_conc[lake_id] = {'best': {ind: round(val, 3) for ind, val in c_best.items()}, 'worst': {ind: round(val, 3) for ind, val in c_worst.items()}}

    ## save data
    with booklet.open(params.lake_acheivable_conc_path, 'n', value_serializer='orjson', key_serializer='uint4', n_buckets=10007) as f:
        for lake_id, val in achieve_conc.items():
            f[lake_id] = val

























































