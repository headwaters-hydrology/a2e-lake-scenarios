#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 12:54:58 2024

@author: mike
"""
import pandas as pd
import hdf5tools
import hdf5plugin
import xarray as xr

import utils, params


pd.options.display.max_columns = 10

##########################################
### Parameters


encoding = {'scale_factor': 1, '_FillValue': -999, 'dtype': 'int16'}


##########################################
### Restructure the input data

def process_periphyton():
    peri_means0 = pd.read_csv(params.peri_mean_csv_path)
    peri_upr0 = pd.read_csv(params.peri_upr_csv_path)

    peri_means1 = peri_means0[[col for col in peri_means0.columns if '_SD' not in col]].copy()
    peri_upr1 = peri_upr0[[col for col in peri_upr0.columns if '_SD' not in col]].copy()

    peri_means1.replace({'Nutrient': {'DIN': 'NO3N'}}, inplace=True)
    peri_upr1.replace({'Nutrient': {'DIN': 'NO3N'}}, inplace=True)

    peri_means2 = peri_means1.set_index(['Nutrient', 'SoF'])
    peri_upr2 = peri_upr1.set_index(['Nutrient', 'SoF', 'UPR'])

    peri_means2['Unshaded_0'] = 0
    peri_means2['Shaded_0'] = 0
    peri_upr2['Unshaded_0'] = 0
    peri_upr2['Shaded_0'] = 0

    peri_means2.columns = pd.MultiIndex.from_tuples([tuple([col.split('_')[0], int(col.split('_')[1])]) for col in peri_means2.columns])
    peri_upr2.columns = pd.MultiIndex.from_tuples([tuple([col.split('_')[0], int(col.split('_')[1])]) for col in peri_upr2.columns])

    peri_means2 = peri_means2.sort_index(axis=1)
    peri_upr2 = peri_upr2.sort_index(axis=1)

    peri_means3 = peri_means2.stack((0, 1))
    peri_means3.name = 'conc'
    peri_means3.index.names = ['indicator', 'source', 'shading', 'biomass']
    peri_means3 = peri_means3.round().astype('int16')
    peri_upr3 = peri_upr2.stack((0, 1))
    peri_upr3.name = 'conc'
    peri_upr3.index.names = ['indicator', 'source', 'upr', 'shading', 'biomass']
    peri_upr3 = peri_upr3.round().astype('int16')

    peri_means4 = peri_means3.to_xarray()
    peri_means4.encoding = encoding
    peri_upr4 = peri_upr3.to_xarray()
    peri_upr4.encoding = encoding

    hdf5tools.xr_to_hdf5(peri_means4.to_dataset(), params.peri_mean_hdf_path)
    hdf5tools.xr_to_hdf5(peri_upr4.to_dataset(), params.peri_upr_hdf_path)































































