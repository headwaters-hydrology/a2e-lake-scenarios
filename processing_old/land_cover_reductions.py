#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 13:11:07 2022

@author: mike
"""
import os
from gistools import vector, rec
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely import intersection
import hdf5tools
import xarray as xr
import dbm
import booklet
import pickle
import zstandard as zstd
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, GradientBoostingRegressor, HistGradientBoostingRegressor, AdaBoostRegressor, ExtraTreesRegressor
from sklearn.preprocessing import OrdinalEncoder
from sklearn.compose import make_column_transformer
from sklearn.compose import make_column_selector
from sklearn.pipeline import make_pipeline

import utils, params

pd.options.display.max_columns = 10

##########################################
### land use/cover

# lcdb_reductions = {'nitrogen': {
#                        'Exotic Forest': 0,
#                        'Forest - Harvested': 0,
#                        'Orchard, Vineyard or Other Perennial Crop': 15,
#                        'Short-rotation Cropland': 30,
#                        'Built-up Area (settlement)': 0,
#                        'High Producing Exotic Grassland': 30,
#                        'Low Producing Grassland': 10,
#                        'Mixed Exotic Shrubland': 0
#                        },
#                    'phosphorus': {
#                        'Exotic Forest': 30,
#                        'Forest - Harvested': 30,
#                        'Orchard, Vineyard or Other Perennial Crop': 50,
#                        'Short-rotation Cropland': 50,
#                        'Built-up Area (settlement)': 0,
#                        'High Producing Exotic Grassland': 30,
#                        'Low Producing Grassland': 10,
#                        'Mixed Exotic Shrubland': 0
#                        },
#                    'sediment': {
#                        'Exotic Forest': 30,
#                        'Forest - Harvested': 30,
#                        'Orchard, Vineyard or Other Perennial Crop': 50,
#                        'Short-rotation Cropland': 50,
#                        'Built-up Area (settlement)': 0,
#                        'High Producing Exotic Grassland': 30,
#                        'Low Producing Grassland': 10,
#                        'Mixed Exotic Shrubland': 0
#                        },
#                    'e.coli': {
#                        'Exotic Forest': 0,
#                        'Forest - Harvested': 0,
#                        'Orchard, Vineyard or Other Perennial Crop': 0,
#                        'Short-rotation Cropland': 30,
#                        'Built-up Area (settlement)': 0,
#                        'High Producing Exotic Grassland': 30,
#                        'Low Producing Grassland': 10,
#                        'Mixed Exotic Shrubland': 0
#                        },
#                    }


# snb_reductions = {'sediment': 30,
#                   'e.coli': 35
#                   }
# dairy_reductions = {'sediment': 65,
#                   'e.coli': 75
#                   }

# snb_reductions = {
#                   'e.coli': 35
#                   }
# dairy_reductions = {
#                   'e.coli': 75
#                   }

features_cols = ['Climate', 'Slope', 'Drainage', 'Wetness']

# param_mapping = {'Visual Clarity': 'sediment',
#                  'E.coli': 'e.coli',
#                  'Dissolved reactive phosporus': 'phosphorus',
#                  'Ammoniacal nitrogen': 'nitrogen',
#                  'Nitrate': 'nitrogen',
#                  'Total nitrogen': 'nitrogen',
#                  'Total phosphorus': 'phosphorus',
#                  'Chlorophyll a': 'e.coli',
#                  'Total Cyanobacteria': 'e.coli',
#                  'Secchi Depth': 'sediment'
#                  }


def land_cover_reductions():
    ### Apply the default reductions
    ## SnB
    # snb_geo = gpd.read_feather(params.snb_geo_clean_path)
    # snb_geo['area_ha'] = snb_geo.geometry.area*0.0001

    # snb_areas = snb_geo.groupby('typology')['area_ha'].sum()

    snb1 = pd.read_csv(params.snb_typo_path)
    snb1['typology'] = snb1.typology.str.title()
    snb1['typology'] = snb1['typology'].str.replace('Bop', 'BoP')

    calc_columns = [col for col in snb1.columns if 'current' in col]
    for col in calc_columns:
        new_col_name = col.replace('load', 'yield')[:-1]+'/ha)'
        snb1[new_col_name] = (snb1[col]/snb1['Area (ha)']).round(3)

    phos2 = ((1 - (snb1['Phosphorus|2035 potential load (kg)'] / snb1['Phosphorus|2015 current load (kg)'])) * 100).round().astype('int8')
    phos2.name = 'phosphorus_reduction'
    nitrate2 = ((1 - (snb1['Nitrogen|2035 potential load (kg)'] / snb1['Nitrogen|2015 current load (kg)'])) * 100).round().astype('int8')
    nitrate2.name = 'nitrogen_reduction'
    sed2 = ((1 - (snb1['Sediment|2035 potential load (t)'] / snb1['Sediment|2015 current load (t)'])) * 100).round().astype('int8')
    sed2.name = 'sediment_reduction'
    snb2 = pd.concat([snb1, phos2, nitrate2, sed2], axis=1).drop('Area (ha)', axis=1)

    # Determine missing typologies for snb
    snb_geo = gpd.read_feather(params.snb_geo_clean_path)

    miss_snb = snb_geo[['typology']][~snb_geo['typology'].isin(snb2.typology)].drop_duplicates(subset=['typology'])

    if not miss_snb.empty:
        raise ValueError('What the heck!')

        snb2['base'] = snb2['typology'].apply(lambda x: x.split('(')[0])
        extra_snb = snb2.groupby('base')[['phosphorus', 'nitrogen']].mean().round().astype('int8')
        snb2.drop('base', axis=1, inplace=True)
        miss_snb['base'] = miss_snb['typology'].apply(lambda x: x.split('(')[0])
        miss_snb1 = pd.merge(miss_snb, extra_snb.reset_index(), on='base').drop('base', axis=1)

        snb2 = pd.concat([snb2, miss_snb1])

    snb3 = snb_geo.merge(snb2, on='typology')

    ## dairy
    dairy1 = pd.read_csv(params.dairy_typo_path)

    calc_columns = [col for col in dairy1.columns if 'current' in col]
    new_col_names = []
    for col in calc_columns:
        new_col_name = col.replace('load', 'yield')[:-1]+'/ha)'
        dairy1[new_col_name] = (dairy1[col]/dairy1['Area (ha)']).round(3)
        new_col_names.append(new_col_name)

    phos2 = ((1 - (dairy1['Phosphorus|2035 potential load (kg)'] / dairy1['Phosphorus|2015 current load (kg)'])) * 100).round().astype('int8')
    phos2.name = 'phosphorus_reduction'
    nitrate2 = ((1 - (dairy1['Nitrogen|2035 potential load (kg)'] / dairy1['Nitrogen|2015 current load (kg)'])) * 100).round().astype('int8')
    nitrate2.name = 'nitrogen_reduction'
    sed2 = ((1 - (dairy1['Sediment|2035 potential load (t)'] / dairy1['Sediment|2015 current load (t)'])) * 100).round().astype('int8')
    sed2.name = 'sediment_reduction'
    dairy2 = pd.concat([dairy1, phos2, nitrate2, sed2], axis=1).drop('Area (ha)', axis=1)
    dairy2 = dairy2.replace({'Wetness': {'Irrig': 'Irrigated'},
                             'Slope': {'Mod': 'Moderate',
                                       'Flat': 'Low'}})
    dairy2['typology'] = dairy2['Climate'] + '/' + dairy2['Slope'] + '/' + dairy2['Drainage'] + '/' + dairy2['Wetness']

    # Determine missing typologies for dairy
    dairy_geo = gpd.read_feather(params.dairy_geo_clean_path)

    d_typo1 = dairy2['typology']

    g_typo1 = pd.Series(dairy_geo['typology'][~dairy_geo['typology'].isin(d_typo1)].unique()).copy()
    t1 = g_typo1.str.split('/')
    typo2 = pd.DataFrame.from_dict(dict(zip(t1.index, t1.values)), orient='index')
    typo2.columns = features_cols
    typo2 = typo2.drop_duplicates()
    typo3 = typo2.copy()

    t2 = pd.Series(dairy_geo['typology'].unique()).str.split('/')
    typo_all = pd.DataFrame.from_dict(dict(zip(t2.index, t2.values)), orient='index')
    typo_all.columns = features_cols

    # The model
    ordinal_encoder = make_column_transformer(
        (
            OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=np.nan),
            make_column_selector(dtype_include="category"),
        ),
        remainder="passthrough",
        # Use short feature names to make it easier to specify the categorical
        # variables in the HistGradientBoostingRegressor in the next step
        # of the pipeline.
        verbose_feature_names_out=False,
    )
    ordinal_encoder.fit(typo_all.astype('category'))

    # model = make_pipeline(
    # ordinal_encoder, HistGradientBoostingRegressor(loss='squared_error', max_iter=100, learning_rate=0.05, early_stopping=False, categorical_features=list(range(4))))
    # model = HistGradientBoostingRegressor(loss='squared_error', categorical_features=list(range(4)))
    # model = RandomForestRegressor()

    train_features = dairy2[features_cols].astype('category').copy()

    est_params = new_col_names.copy()
    est_params.extend(['phosphorus_reduction', 'nitrogen_reduction', 'sediment_reduction'])

    for param in est_params:
        model = GradientBoostingRegressor(loss='absolute_error')

        train_targets = dairy2[param].copy()

        model.fit(ordinal_encoder.transform(train_features), train_targets.values)

        typo3[param] = model.predict(ordinal_encoder.transform(typo2.astype('category'))).round(2)

    # Export the results for checking
    typo3.to_csv(params.dairy_model_typo_path, index=False)

    # Combine with original data
    dairy3 = pd.concat([dairy2, typo3])
    dairy3['typology'] = dairy3['Climate'] + '/' + dairy3['Slope'] + '/' + dairy3['Drainage'] + '/' + dairy3['Wetness']

    dairy4 = dairy_geo.merge(dairy3.drop(features_cols, axis=1), on='typology')

    drop_cols = [col for col in dairy4.columns if 'load' in col]
    combo1 = pd.concat([snb3, dairy4]).reset_index(drop=True).drop(drop_cols, axis=1).drop(['Sediment|2015 current yield (t/ha)', 'sediment_reduction'], axis=1)
    combo1 = combo1.rename(columns={'Phosphorus|2015 current yield (kg/ha)': 'phosphorus_yield', 'Nitrogen|2015 current yield (kg/ha)': 'nitrogen_yield'})

    # combo2 = combo1.drop(set(param_mapping.values()), axis=1)
    # combo2 = combo1.rename(columns={'nitrogen': 'total nitrogen', 'phosphorus': 'total phosphorus'})

    utils.gpd_to_feather(combo1, params.snb_dairy_red_path)

    combo3 = combo1.drop(['geometry'], axis=1).drop_duplicates(subset=['typology'])
    lcdb_extras = combo3.drop(['farm_type', 'typology'], axis=1).groupby('land_cover').mean().round(2)

    ### LCDB - taken care of in other script
    lcdb0 = gpd.read_feather(params.lcdb_clean_path)

    ## Reductions
    lcdb_red = pd.read_csv(params.lcdb_reductions_csv_path).drop(['sediment_reduction', 'e.coli_reduction'], axis=1)

    lcdb1 = lcdb0.merge(lcdb_red, on='typology', how='outer').reset_index(drop=True)
    lcdb1.loc[lcdb1['nitrogen_reduction'].isnull(), ['nitrogen_reduction', 'phosphorus_reduction']] = 0

    ## Yields
    lcdb_yields = pd.read_csv(params.lcdb_yields_csv_path).drop(['sediment_yield'], axis=1)

    lcdb2 = lcdb1.merge(lcdb_yields, on='typology', how='outer').reset_index(drop=True)
    lcdb2.loc[lcdb2['nitrogen_yield'].isnull(), ['nitrogen_yield', 'phosphorus_yield']] = 0

    lcdb2 = lcdb2.dropna() # Remove the native vegetation class - adding later...

    lcdb2.loc[lcdb2.land_cover == 'Low Producing Grassland', lcdb_extras.columns] = lcdb_extras.loc['Sheep and Beef'].values
    lcdb2.loc[lcdb2.land_cover == 'High Producing Exotic Grassland', lcdb_extras.columns] = lcdb_extras.loc['Dairy'].values

    utils.gpd_to_feather(lcdb2, params.lcdb_red_path)

    ## Save the reductions only
    lcdb3 = lcdb2.drop('geometry', axis=1).drop_duplicates(subset=['typology'])

    lc_red = pd.concat([combo3, lcdb3])
    lc_red['new_categories'] = lc_red['land_cover'].replace(params.combine_land_covers)

    lc_red.to_csv(params.lc_red_csv_path, index=False)






































































































