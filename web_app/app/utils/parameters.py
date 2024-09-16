#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 08:38:28 2023

@author: mike
"""
import os
import pathlib
import numpy as np

#################################################
#### Global parameters

### Paths
assets_path = pathlib.Path(os.path.realpath(os.path.dirname(__file__))).parent.joinpath('assets')

app_base_path = pathlib.Path('/assets')

base_data_url = 'https://b2.tethys-ts.xyz/file/'

# lc_url = '{}olw-data/olw-sc008/olw_land_cover_reductions.gpkg'.format(base_data_url)
# lakes_red_url = '{}olw-data/olw-sc008/olw_lakes_reductions.csv.zip'.format(base_data_url)

lakes_points_path = app_base_path.joinpath('lake_points.pbf')
lakes_poly_path = assets_path.joinpath('lake_poly_gbuf.blt')
lakes_catch_reaches_path = assets_path.joinpath('lake_reaches_gbuf.blt')
# lakes_sites_names_path = assets_path.joinpath('lake_wq_sites_names.blt')
lakes_catch_path = assets_path.joinpath('lake_catch_major_gbuf.blt')

# lakes_sites_catch_reductions_path = assets_path.joinpath('sites_area_reductions.blt')
lakes_lc_red_yields_path = assets_path.joinpath('lake_land_cover_reductions_yields.blt')

lakes_ref_conc_path = assets_path.joinpath('lake_reference_conc.blt')

lakes_model_medians_path = assets_path.joinpath('lake_model_conc_median.blt')
lakes_moni_medians_path = assets_path.joinpath('lake_moni_conc_median.blt')
lakes_moni_data_path = assets_path.joinpath('lake_moni_conc_data.blt')

## misc params
lu_order = ['Dairy', 'Sheep and Beef', 'Short-rotation Crop', 'Perennial Crop', 'Forestry', 'Native Vegetation', 'Urban', 'Other']

# ecoli_reductions = {
#     'Sheep and Beef': 0.5,
#     'Dairy': 0.5,
#     'Other': 0,
#     'Forestry': 0,
#     'Short-rotation Crop': 0,
#     'Perennial Crop': 0,
#     'Urban': 0,
#     'Native Vegetation': 0
#     }

lakes_indicator_dict = {
    'TN': 'Total nitrogen',
    'TP': 'Total phosphorus',
    'CHLA': 'Chlorophyll a',
    'Secchi': 'Secchi depth'
    }

indicator_rounding = {
    'TN': 0,
    'TP': 0,
    'CHLA': 1,
    'Secchi': 2
    }

indicator_str_format = {
    'TN': '{:.0f}',
    'TP': '{:.0f}',
    'CHLA': '{:.1f}',
    'Secchi': '{:.2f}'
    }

lakes_reduction_cols = list(lakes_indicator_dict.values())

nps_mapping = {
    'TN': 'Total nitrogen',
    'TP': 'Total phosphorus',
    'CHLA': 'Chla',
    # 'ECOLI': 'E.coli',
    }

tbl_cols_dict = {
    'land_cover': ['', 'Land use'],
    'land_area': ['Current', 'Land area %'],
    'load_n': ['Current', 'Nitrogen Load %'],
    'load_p': ['Current', 'Phosphorus Load %'],
    'nitrogen': ['Scenario', 'Nitrogen Mitigation %'],
    'phosphorus': ['Scenario', 'Phosphorus Mitigation %'],
    'new_land_area': ['Scenario', 'Land area %'],
    # 'new_load': ['Scenario', 'Load contribution %']
    }

ind_tbl_cols_main_dict = {
    'name': 'Name',
    'conc': 'Median concentration (mg/m³)',
    }

ind_tbl_cols_secchi_dict = {
    'name': 'Name',
    'conc': 'Median depth (m)',
    }

fill_opacity = 0.5
plot_colors = {
    'Band A': f'rgba(204, 235, 197, {fill_opacity})',
    'Band B': f'rgba(255, 255, 204, {fill_opacity})',
    'Band C': f'rgba(254, 217, 166, {fill_opacity})',
    'Band D': f'rgba(251, 180, 174, {fill_opacity})',
    'Current': 'rgb(55, 126, 184)',
    'Scenario': 'rgb(166, 86, 40)',
    }

### Layout
space_between_content = 10
header_font_size = 17
map_height = '47vh' # for combined website
center = [-41.1157, 172.4759]
zoom = 6

hovercard_width = 300
hovercard_open_delay = 1000

attribution = 'Map data &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'













