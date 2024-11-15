#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 08:12:43 2024

@author: mike
"""
import os
import sys
import pathlib
import geopandas as gpd
import pandas as pd
import numpy as np
import booklet
from shapely.ops import unary_union
from shapely import intersection, difference, intersects, symmetric_difference, make_valid

pd.options.display.max_columns = 10

#########################################################
### Parameters

script_path = pathlib.Path(os.path.realpath(os.path.dirname(__file__)))

processing_path = script_path.parent.joinpath('processing')

sys.path.append(str(processing_path))

import params, utils

orig_lu_path = script_path.joinpath('original_land_use.gpkg')
# wise_lu_path = script_path.joinpath('wise_land_use.gpkg')
# combo_lu_path = script_path.joinpath('combo_land_use.gpkg')

orig_lu_reductions_path = script_path.joinpath('lake_reductions_41314.gpkg')

# site_id = 'ew-00108'
# nzsegment = 3047194

wise_lu_mapping = {
    'Area outside region': 'Other',
    'Commercial': 'Urban',
    'Community Services': 'Urban',
    'Freshwater': 'Other',
    'Manufacturing': 'Urban',
    'Mines and Quarries': 'Other',
    'Utilities': 'Urban',
    'Vacant - urban land': 'Urban',
    'Dairying': 'Dairy',
    'Forestry': 'Forestry',
    'Horticulture': 'Perennial Crop',
    'Indigenous': 'Native Vegetation',
    'Lifestyle': 'Other',
    'Low Density Residential': 'Urban',
    'Medium to High Density Residential': 'Urban',
    'Other Agriculture': 'Sheep and Beef',
    'Other Cropping': 'Short-rotation Crop',
    'Other Exotic': 'Forestry',
    'Sheep and Beef': 'Sheep and Beef',
    'Urban Parks and Recreation': 'Urban',
    'Vegetable Cropping': 'Vegetable Crop',
    'Wetlands': 'Native Vegetation',
    }

######################################################
### Clip out the wise land use


# with booklet.open(params.sites_catch_major_path) as f:
#     catch = f[nzsegment]

# poly = gpd.GeoSeries(data=[catch], crs=4326).to_crs(2193)[0]
# wise0 = gpd.read_feather(params.wise_lu_cleaned_path)

# wise1 = wise0[wise0.geometry.intersects(poly)].copy()

# geo = intersection(wise1.geometry.tolist(), poly)

# wise1['geometry'] = geo
# wise2 = wise1.loc[~wise1.geometry.is_empty, ['LU_Name2023', 'geometry']].rename(columns={'LU_Name2023': 'land_use'}).copy()

# wise3 = wise2.replace({'land_use': wise_lu_mapping}).dissolve('land_use').reset_index()

# wise3.to_file(wise_lu_path)


####################################################
### Save other land use layers

# with booklet.open(params.sites_land_cover_reductions_path) as f:
#     lu_catch = f[nzsegment][['land_cover', 'geometry']].rename(columns={'land_cover': 'land_use'}).copy()

# geo = intersection(lu_catch.geometry.tolist(), poly)

# lu_catch['geometry'] = geo
# lu_catch1 = lu_catch.loc[~lu_catch.geometry.is_empty].copy()

# lu_catch2 = lu_catch1.replace({'land_use': params.combine_land_covers}).dissolve('land_use').reset_index()

# lu_catch2.to_file(combo_lu_path)



lu_catch = gpd.read_file(orig_lu_reductions_path)[['land_cover', 'geometry']].rename(columns={'land_cover': 'land_use'}).copy()

lu_catch2 = lu_catch.replace({'land_use': params.combine_land_covers}).dissolve('land_use').reset_index()

lu_catch2.to_file(orig_lu_path)



















































































