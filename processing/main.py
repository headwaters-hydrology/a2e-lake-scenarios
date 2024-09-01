#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 12:57:57 2024

@author: mike
"""
from lakes_data_monitored import lakes_monitored_conc
from lakes_model_conc import process_model_conc_lakes
from lakes_catch_delineation import delin_lakes_catch_rec
from lakes_geo_processing import lakes_points_poly_process
from lakes_inflow_loads import process_inflow_loads_lakes
from lakes_reference_conc import process_reference_conc_rec
from lcdb_processing import lcdb_processing
from land_cover_reductions import land_cover_reductions
from combine_yields_reductions_for_app import combine_red_yields

####################################################
### Order of processing operations

## Lake monitored WQ data
lakes_monitored_conc()

## Lakes modelled conc
process_model_conc_lakes()

## Lake catch delineation
delin_lakes_catch_rec()

## Lakes points and polygons
lakes_points_poly_process()

## Lake inflows
process_inflow_loads_lakes()

## Estimate reference conc at all relevant lake segments
process_reference_conc_rec()

## LCDB processing
lcdb_processing()

## land cover reductions nationally
land_cover_reductions()

## Land cover reductions per site/catchment - cmd
# lakes_land_cover_reductions.py

## Split the land cover yields to the catchments above sites - cmd

# LCDB
# lakes_land_cover_lcdb_yields.py

# Farmland
# lakes_land_cover_farm_yields.py

# Combine
# lakes_land_cover_combine_yields.py

## Combine yields and reductions
combine_red_yields()

## Make ts data stats for app
wq_data_stats()

## Make site catchment names for app
sites_assign_catch_names()

## Process perphyton data
process_periphyton()

## Process ecoli yields
ecoli_yields()


















































