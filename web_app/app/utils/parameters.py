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
# rivers_red_url = '{}olw-data/olw-sc008/olw_rivers_reductions.csv.zip'.format(base_data_url)

## Rivers montiroing sites scenarios
rivers_sites_path = app_base_path.joinpath('river_wq_sites.pbf')
rivers_sites_ind_path = assets_path.joinpath('river_wq_sites_ind.blt')
rivers_sites_reaches_path = assets_path.joinpath('sites_reaches_gbuf.blt')
rivers_sites_names_path = assets_path.joinpath('river_wq_sites_names.blt')
rivers_sites_catch_path = assets_path.joinpath('sites_catch_major_gbuf.blt')

# rivers_sites_catch_reductions_path = assets_path.joinpath('sites_area_reductions.blt')
sites_lc_red_yields_path = assets_path.joinpath('sites_land_cover_reductions_yields.blt')

rivers_sites_ref_conc_path = assets_path.joinpath('sites_reference_conc.blt')

rivers_sites_wq_stats_path = assets_path.joinpath('river_wq_data_stats.blt')
rivers_data_app_path = assets_path.joinpath('river_wq_data_for_app.blt')

sites_catch_names_path = assets_path.joinpath('sites_catch_names.blt')

## misc params
lu_order = ['Dairy', 'Sheep and Beef', 'Short-rotation Crop', 'Perennial Crop', 'Forestry', 'Native Vegetation', 'Urban', 'Other']

ecoli_reductions = {
    'Sheep and Beef': 0.5,
    'Dairy': 0.5,
    'Other': 0,
    'Forestry': 0,
    'Short-rotation Crop': 0,
    'Perennial Crop': 0,
    'Urban': 0,
    'Native Vegetation': 0
    }

rivers_indicator_dict = {
    'ECOLI': 'E. coli',
    'DRP': 'Dissolved reactive phosphorus',
    'NO3N': 'Nitrate nitrogen',
    'TN': 'Total nitrogen',
    'TP': 'Total phosphorus'
    }

rivers_reduction_cols = list(rivers_indicator_dict.values())

nps_mapping = {
    'DRP': 'DRP',
    # 'ECOLI': 'E.coli',
    'NO3N': 'Nitrate'
    }

tbl_cols_dict = {
    'land_cover': ['', 'Land use'],
    'land_area': ['Current', 'Land area %'],
    'load': ['Current', 'Load contribution %'],
    'mitigation': ['Scenario', 'Mitigation effectiveness %'],
    'new_land_area': ['Scenario', 'Land area %'],
    'new_load': ['Scenario', 'Load contribution %']
    }

ind_tbl_cols_main_dict = {
    'name': 'Name',
    'conc': 'Median concentration (mg/l)',
    }

ind_tbl_cols_ecoli_dict = {
    'name': 'Name',
    'conc': 'Median concentration (CFU/100ml)',
    }

peri_tbl_cols_dict = {
    'name': ['', 'Name'],
    'conc_shaded': ['Best estimate Q92 Chla (mg/m²)', 'Shaded'],
    'conc_unshaded': ['Best estimate Q92 Chla (mg/m²)', 'Unshaded'],
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

upr = [5, 10, 15, 20, 25, 30, 50]

plot_colors_upr = {
    'Band A': f'rgba(204, 235, 197, {fill_opacity})',
    'Band B': f'rgba(255, 255, 204, {fill_opacity})',
    'Band C': f'rgba(254, 217, 166, {fill_opacity})',
    'Band D': f'rgba(251, 180, 174, {fill_opacity})',
    'Scenario Shaded': 'rgb(166,97,26)',
    'Scenario Unshaded': 'rgb(223,194,125)',
    'Current Shaded': 'rgb(5,113,176)',
    'Current Unshaded': 'rgb(146,197,222)',
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













