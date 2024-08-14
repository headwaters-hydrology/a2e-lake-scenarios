#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 12:57:57 2024

@author: mike
"""
from lawa_wq_sites_for_app import lawa_wq_sites_for_app
from site_catch_delineation import delin_sites_catch_rec
# from rivers_process_flows import process_flows_rec
# from rivers_process_loads import process_loads_rec
from site_reference_conc import process_reference_conc_rec
#  catch_land_cover_lcdb_reductions
# import land_cover_per_rec_catch
from lcdb_processing import lcdb_processing
from land_cover_reductions import land_cover_reductions
from sites_stats import wq_data_stats
from combine_yields_reductions_for_app import combine_red_yields
from sites_catch_names import sites_assign_catch_names
from periphyton_data import process_periphyton
from ecoli_spatial_data_yields import ecoli_yields


####################################################
### Order of processing operations - old...

## Site catch delineation
# delin_sites_catch_rec()

## Estimate flow differences at all relevant segments
# process_flows_rec()

## Estimate loads at all relevant segments and sites
# process_loads_rec()

## Estimate reference conc at all relevant segments and sites
# process_reference_conc_rec()

## Split all of the land cover/use reductions polygons to the site catchments
# catch_land_cover_reductions()

## Estimate the land cover loads from the river loads (multiprocessing)
# land_cover_per_rec_catch.py


####################################################
### Order of processing operations - new...

## Site WQ data
lawa_wq_sites_for_app()

## Site catch delineation
delin_sites_catch_rec()

## Estimate reference conc at all relevant segments and sites
process_reference_conc_rec()

## LCDB processing
lcdb_processing()

## land cover reductions nationally
land_cover_reductions()

## Land cover reductions per site/catchment - cmd
# catch_land_cover_reductions.py

## Split the land cover yields to the catchments above sites - cmd

# LCDB
# catch_land_cover_lcdb_yields.py

# Farmland
# catch_land_cover_farm_yields.py

# Combine
# catch_land_cover_combine_yields.py

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


















































