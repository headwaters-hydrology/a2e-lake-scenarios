#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 00:42:44 2024

@author: mike
"""
import os
import pathlib


###############################################
### Global parameters

script_path = pathlib.Path(os.path.realpath(os.path.dirname(__file__)))
data_path = script_path.parent.joinpath('data')

if not data_path.is_dir():
    raise ValueError(f'{data_path} does not exist')

nzrec_data_path = '/home/mike/git/nzrec/data'
rec_rivers_feather = data_path.joinpath('rec25_rivers_clean.feather')

lakes_poly_path = data_path.joinpath('lake_polygons_fenz.gpkg')
lakes_catch_path = data_path.joinpath('FENZ_v1_LakeCatchments.shp')

rc_poly_path = data_path.joinpath('regional-council-2023-generalised.gpkg')

### Loads/concs
indicators = ['CHLA', 'NH4N', 'Secchi', 'TN', 'TP', 'ECOLI']

lake_inflows_gpkg = data_path.joinpath('lake_inflows_segments.gpkg')
lake_outflows_gpkg = data_path.joinpath('lake_outflows_segments.gpkg')
lake_outflows_flows_csv = data_path.joinpath('lake_outflows_flows.csv')
lake_catch_segs_gpkg = data_path.joinpath('lake_catchment_segments.gpkg')

lake_inflows_blt = data_path.joinpath('lake_inflows_segments.blt')
lake_outflows_blt = data_path.joinpath('lake_outflows_segments.blt')
lake_catch_segs_blt = data_path.joinpath('lake_catchment_segments.blt')

lake_moni_data_raw = data_path.joinpath('lawa_lake_wq_field_samples_ts.feather')
lake_site_data_raw = data_path.joinpath('lawa_lake_wq_field_samples_sites.gpkg')
lake_mtypes_data_raw = data_path.joinpath('lawa_lake_wq_field_samples_sites_mtypes.csv')

lake_moni_data = data_path.joinpath('lawa_lake_wq_field_samples_ts_filtered.feather')
lake_site_data = data_path.joinpath('lawa_lake_wq_field_samples_sites_filtered.gpkg')

rivers_loads_csv = data_path.joinpath('PredictedYeildDF_REC2_20July2023.csv')

lake_load_inflows = data_path.joinpath('lake_load_inflows.blt')
lake_model_conc_csv = data_path.joinpath('LakeRF_WQModel_Predictions04Feb2022.csv')
lake_model_conc_blt = data_path.joinpath('lake_model_conc_median.blt')
lake_moni_conc_blt = data_path.joinpath('lake_moni_conc_median.blt')
lake_moni_conc_csv = data_path.joinpath('lake_moni_conc_median.csv')
lake_moni_conc_data_blt = data_path.joinpath('lake_moni_conc_data.blt')
lake_metadata_blt = data_path.joinpath('lake_metadata.blt')

lake_lc_data_csv = data_path.joinpath('lake_land_cover_areas_yields.csv')
lake_lc_summ_csv = data_path.joinpath('lake_land_cover_summary.csv')

# river_flows_rec_path = data_path.joinpath('NZRiverMaps_hydrology_2023-01-09.csv')
# rivers_conc_csv_path1 = data_path.joinpath('NutrientConcsYields.csv')
# rivers_conc_csv_path2 = data_path.joinpath('EcoliConcsYields.csv')
# rivers_conc_csv_path3 = data_path.joinpath('updated-suspended-sediment-yield-estimator-and-estuarine-tra.csv')

# river_flows_rec_diff_path = data_path.joinpath('rec_median_flows_diff.blt')
# river_loads_rec_diff_path = data_path.joinpath('rec_loads_diff.blt')
# river_loads_rec_diff_csv_path = data_path.joinpath('rec_loads_diff.csv.zip')
# site_catch_rec_loads_path = data_path.joinpath('sites_catch_loads.blt')

## Sites, indicators and ts data
# wq_sites_raw_path = data_path.joinpath('lawa_river_wq_field_samples_sites.gpkg')
# wq_data_raw_path = data_path.joinpath('lawa_river_wq_field_samples_ts.feather')

# sites_remove_path = data_path.joinpath('removed_sites.csv')

# wq_sites_gbuf_path = data_path.joinpath('river_wq_sites.pbf')
# wq_sites_gpkg_path = data_path.joinpath('river_wq_sites.gpkg')
# wq_sites_ind_path = data_path.joinpath('river_wq_sites_ind.blt')
# wq_data_path = data_path.joinpath('river_wq_data.blt')
# wq_data_app_path = data_path.joinpath('river_wq_data_for_app.blt')
# wq_sites_names_path = data_path.joinpath('river_wq_sites_names.blt')

## Catchment delineation
lakes_reach_mapping_path = data_path.joinpath('lake_reaches_mapping.blt')
lakes_reach_gbuf_path = data_path.joinpath('lake_reaches_gbuf.blt')
lakes_catch_major_path = data_path.joinpath('lake_catch_major.blt')
lakes_catch_major_gbuf_path = data_path.joinpath('lake_catch_major_gbuf.blt')

lakes_catch_gpkg_path = data_path.joinpath('lake_catch_major.gpkg')
lakes_catch_minor_path = data_path.joinpath('lake_catch_minor.blt')

lakes_catch_names_path = data_path.joinpath('lake_catch_names.blt')

limnotrack_path = data_path.joinpath('limnotrack')

lt_reach_data_csv = limnotrack_path.joinpath('reach_data.csv')
lt_inflows_gpkg = limnotrack_path.joinpath('lake_inflows_limnotrack.gpkg')
lt_outflows_gpkg = limnotrack_path.joinpath('lake_outflows_limnotrack.gpkg')

## Lakes polygons and points
lakes_poly_gbuf_path = data_path.joinpath('lake_poly_gbuf.blt')
lakes_points_gbuf_path = data_path.joinpath('lake_points.pbf')

lakes_poly_gpkg_path = data_path.joinpath('lake_poly.gpkg')
lakes_points_gpkg_path = data_path.joinpath('lake_points.gpkg')

## Reference conc
lake_ref_conc_csv = data_path.joinpath('LakeReferenceIndicators.csv')
lake_missing_ref_csv = data_path.joinpath('lake_ref_missing.csv')
rivers_ref_conc3_csv_path = data_path.joinpath('reference_conc_rec_level_3.csv')
rivers_ref_conc2_csv_path = data_path.joinpath('reference_conc_rec_level_2.csv')
rivers_ref_conc_csv_path = data_path.joinpath('reference_conc_rec_clean.csv.zip')
rec_reference_conc_path = data_path.joinpath('rec_reference_conc.blt')
lake_reference_conc_path = data_path.joinpath('lake_reference_conc.blt')

## NPS-FM stats and bands
wq_data_stats_path = data_path.joinpath('river_wq_data_stats.blt')

nps_mapping = {
    'DRP': 'DRP',
    # 'ECOLI': 'E.coli',
    'NO3N': 'Nitrate'
    }

### Land cover
lu_order = ['Dairy', 'Sheep and Beef', 'Short-rotation Crop', 'Perennial Crop', 'Forestry', 'Native Vegetation', 'Urban', 'Other']

combine_land_covers = {
    'Low Producing Grassland': 'Sheep and Beef',
    'High Producing Exotic Grassland': 'Dairy',
    'Orchard, Vineyard or Other Perennial Crop': 'Perennial Crop',
    'Short-rotation Cropland': 'Short-rotation Crop',
    'Exotic Forest': 'Forestry',
    'Mixed Exotic Shrubland': 'Sheep and Beef',
    'Forest - Harvested': 'Forestry',
    'Built-up Area (settlement)': 'Urban',
    'Estuarine Open Water': 'Other',
    'Gravel or Rock': 'Other',
    'Lake or Pond': 'Other',
    # 'Mangrove': 'Other',
    'Not land': 'Other',
    'River': 'Other',
    'Sand or Gravel': 'Other',
    'Urban Parkland/Open Space': 'Other',
    'Surface Mine or Dump': 'Other',
    'Transport Infrastructure': 'Other',
    'Permanent Snow and Ice': 'Other',
    'Herbaceous Saline Vegetation': 'Other',
    'Herbaceous Freshwater Vegetation': 'Other',
    'Depleted Grassland': 'Sheep and Beef',
    }


# typo_red_corrections = {
#     3076139: {
#         'Warm/Low/Well/Moist': {'phosphorus_reduction': 36, 'nitrogen_reduction': 39},
#         'Cool/Low/Well/Moist': {'phosphorus_reduction': 34, 'nitrogen_reduction': 32},
#         'Cool/Low/Light/Moist': {'phosphorus_reduction': 35, 'nitrogen_reduction': 35},
#         'Warm/Low/Light/Moist': {'phosphorus_reduction': 35, 'nitrogen_reduction': 35},
#         'Warm/Moderate/Light/Moist': {'phosphorus_reduction': 35, 'nitrogen_reduction': 35},
#         'Warm/Moderate/Well/Moist': {'phosphorus_reduction': 35, 'nitrogen_reduction': 35},
#         'High Producing Exotic Grassland': {'phosphorus_reduction': 35, 'nitrogen_reduction': 35},
#         }
#     }

yields_csv_path = data_path.joinpath('LandUseToLossRatesV0.csv')
yields_gpkg_path = data_path.joinpath('SrinivasanTypes.gpkg')

lcdb_reductions_csv_path = data_path.joinpath('lcdb_reductions.csv')
lcdb_yields_csv_path = data_path.joinpath('lcdb_yields.csv')

snb_geo_path = data_path.joinpath('SnB_Typologies.shp')
dairy_geo_path = data_path.joinpath('Dairy_Typologies.shp')
dairy_geo_clean_path = data_path.joinpath('dairy_typologies_clean.feather')
snb_geo_clean_path = data_path.joinpath('snb_typologies_clean.feather')

snb_typo_path = data_path.joinpath('typologies to reductions - snb.csv')
dairy_typo_path = data_path.joinpath('typologies to reductions - dairy.csv')
dairy_model_typo_path = data_path.joinpath('dairy_modelled_typologies.csv')

snb_dairy_red_path = data_path.joinpath('snb_dairy_reductions.feather')

lc_red_csv_path = data_path.joinpath('typology_reductions.csv')

lcdb_path = data_path.joinpath('lcdb_v5.fgb')
lcdb_clean_path = data_path.joinpath('lcdb_clean.feather')
lcdb_red_path = data_path.joinpath('lcdb_reductions.feather')
lcdb_yields_path = data_path.joinpath('lcdb_yields.feather')
snb_dairy_red_path = data_path.joinpath('snb_dairy_reductions.feather')

awm_red_path = data_path.joinpath('awm_reductions.feather')
awm_yields_path = data_path.joinpath('awm_yields.feather')

lakes_farm_yields_blt_path = data_path.joinpath('lake_land_cover_farm_yields.blt')
lakes_lcdb_yields_blt_path = data_path.joinpath('lake_land_cover_lcdb_yields.blt')
lakes_yields_blt_path = data_path.joinpath('lake_land_cover_yields.blt')

lakes_land_cover_reductions_path = data_path.joinpath('lake_land_cover_reductions.blt')
lakes_lc_red_yields_path = data_path.joinpath('lake_land_cover_reductions_yields.blt')
lakes_lc_red_yields_gpkg = data_path.joinpath('lake_land_cover_reductions_yields.gpkg')


lakes_land_cover_catch_path = data_path.joinpath('lake_land_cover_per_catch.blt')
example_catch_land_cover_path = data_path.joinpath('example_land_cover_base.gpkg')
example_catch_land_cover_rec_path = data_path.joinpath('example_land_cover_by_rec_watershed.gpkg')
sites_land_cover_catch_loads_path = data_path.joinpath('sites_land_cover_per_catch_loads.blt')
sites_land_cover_catch_loads_feather_path = data_path.joinpath('sites_land_cover_per_catch_loads.feather')
sites_land_cover_catch_loads_gpkg_path = data_path.joinpath('sites_land_cover_per_catch_loads.gpkg')

land_cover_loads_agg_csv_path = data_path.joinpath('land_cover_yields_agg.csv')

sites_lc_csv = data_path.joinpath('sites_land_cover_areas_yields.csv')

sites_snb_dairy_yields_gpkg_path = data_path.joinpath('sites_snb_dairy_yields.gpkg')
sites_snb_dairy_yields_csv_path = data_path.joinpath('sites_snb_dairy_yields.csv')

lc_summ_csv_path = data_path.joinpath('land_cover_summary.csv')

### Ton model TN data comparisons
ton_tn_loads_csv_path = data_path.joinpath('ton_tn_yield_model.csv')

ton_lc_yields = {'Bare': 0,
                 'Water': 0,
                 'Cropland': 11,
                 'Dairy_Dry': 28,
                 'Dairy_Irrigated': 30,
                 'Dairy_Moist': 17,
                 'Dairy_Wet': 37,
                 'Forestry': 8,
                 'Native': 2,
                 'OrchardVineyard': 19,
                 'SheepBeef_Flat': 8,
                 'SheepBeef_Hill': 4,
                 'Urban': 10
                 }

# ton_lc_map = {'Bare': 'Unchangeable',
#               'Cropland': 'Cropland',
#               'Dairy_Dry': 'Dairy',
#               'Dairy_Irrigated': 'Dairy',
#               'Dairy_Moist': 'Dairy',
#               'Dairy_Wet': 'Dairy',
#               'Forestry': 'Exotic Forest/Shrubland',
#               'Native': 'Unimprovable',
#               'OrchardVineyard': 'Cropland',
#               'SheepBeef_Flat': 'Sheep and Beef',
#               'SheepBeef_Hill': 'Sheep and Beef',
#               'Urban': 'Unimprovable',
#               'Water': 'Unchangeable',
#               }


### Periphyton
peri_path = data_path.joinpath('periphyton')
peri_mean_csv_path = peri_path.joinpath('MeanValue_AllNutrientShadeCriteria.csv')
peri_upr_csv_path = peri_path.joinpath('AllNutrientCriteriaGLMNumeric_25Sept2023.csv')

peri_mean_hdf_path = peri_path.joinpath('periphyton_means.h5')
peri_upr_hdf_path = peri_path.joinpath('periphyton_upr.h5')


### E.coli
ecoli_reductions = {
    'Sheep and Beef': 0.5,
    'Dairy': 0.5,
    'Other': 0,
    'Exotic Forest': 0,
    'Short-rotation Crop': 0,
    'Perennial Crop': 0,
    'Urban': 0,
    'Native Vegetation': 0
    }

# lu_to_land_cover = {
#     'SheepBeef': 'Sheep and Beef',
#     'Dairy': 'Dairy',
#     'OrchardVineyard': 'Cropland',
#     'Cropland': 'Cropland',
#     'Forestry': 'Exotic Forest',
#     'Urban': 'Urban',
#     'Bare': 'Other',
#     'Water': 'Other',
#     'Native': 'Native Vegetation',
#     }


ecoli_path = data_path.joinpath('ecoli')
ecoli_drain_tif_path = ecoli_path.joinpath('Drainage.tif')
ecoli_elev_tif_path = ecoli_path.joinpath('Elevation.tif')
ecoli_lu_tif_path = ecoli_path.joinpath('LandUse.tif')
ecoli_lu_csv_path = ecoli_path.joinpath('LandUseLookup.csv')
ecoli_input_params = ecoli_path.joinpath('DetailedParametersECOLIconcentrationModel_8.csv')

ecoli_catch_yields_path = ecoli_path.joinpath('ecoli_catch_yields.blt')






