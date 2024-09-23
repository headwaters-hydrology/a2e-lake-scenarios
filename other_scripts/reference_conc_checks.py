#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 18:10:58 2024

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
### Parameters

indicators = ('TN', 'TP', 'CHLA', 'Secchi')


#################################################
### Checks

above = {ind: 0 for ind in indicators}
below = {ind: 0 for ind in indicators}
between = {ind: 0 for ind in indicators}
total = 0

res_list = []
with booklet.open(params.lake_metadata_blt) as meta:
    with booklet.open(params.lake_model_moni_conc_blt) as current:
        with booklet.open(params.lake_acheivable_conc_path) as scenario:
            with booklet.open(params.lake_reference_conc_path) as ref:
                for lake_id in meta:
                    if (lake_id in scenario) and (lake_id in ref) and (lake_id in current):
                        current0 = current[lake_id]
                        scenario0 = scenario[lake_id]
                        ref0 = ref[lake_id]

                        for ind in indicators:
                            c1 = current0[ind]
                            c1_diff = c1*0.2

                            s1 = scenario0[ind]
                            ref1 = ref0[ind]

                            if (ind != 'Secchi') and (ref1 > c1) and ((c1 - c1_diff) > s1):
                                above[ind] += 1
                            elif (ind == 'Secchi') and (ref1 < c1) and ((c1 + c1_diff) < s1):
                                above[ind] += 1

                            if (ind != 'Secchi') and (ref1 < (s1 - c1_diff)):
                                below[ind] += 1
                            elif (ind == 'Secchi') and (ref1 > (s1 + c1_diff)):
                                below[ind] += 1

                            res_list.append([lake_id, ind, c1, s1, ref1])

                        total += 1


above_perc = {}
for ind, a in above.items():
    above_perc[ind] = round((a/total) * 100)

below_perc = {}
for ind, a in below.items():
    below_perc[ind] = round((a/total) * 100)


res = pd.DataFrame(res_list, columns=['LFENZID', 'indicator', 'current', 'scenario', 'reference'])

res2 = res.set_index(['LFENZID', 'indicator']).round(2)

res2.to_csv('/home/mike/git/a2e-lake-scenarios/data/lakes_results_comparisons.csv')














































