#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 14:27:59 2024

@author: mike
"""
import booklet
import pandas as pd
import utils, params





#################################################
### Parameters


catch_id = 15315077
pres_path = params.data_path.joinpath('presentation')

land_use_mataura_gpkg = pres_path.joinpath('land_use_mataura.gpkg')

################################################
### Processing

awm_red = pd.read_feather(params.awm_red_path).set_index('land_cover')
awm_yields = pd.read_feather(params.awm_yields_path).set_index('land_cover')

awm = pd.concat([awm_yields, awm_red], axis=1)
awm['area_ha'] = 0

combo_dict = {}
with booklet.open(params.sites_land_cover_reductions_path) as r:
    with booklet.open(params.sites_yields_blt_path) as y:
        red = r[catch_id]
        yields = y[catch_id]
        if 'farm_type' in yields.columns:
            yields = yields.drop('farm_type', axis=1)

        combo0 = red.merge(yields.drop(['geometry', 'area_ha'], axis=1), on=['land_cover', 'typology']).drop('typology', axis=1)
        combo1 = combo0.replace({'land_cover': params.combine_land_covers})
        combo2 = combo1[['land_cover', 'area_ha', 'geometry']].dissolve('land_cover', aggfunc='sum')
        values = utils.area_weighted_mean(combo1.reset_index(), 'land_cover', ['phosphorus_yield', 'nitrogen_yield', 'phosphorus_reduction', 'nitrogen_reduction']).round(3)
        combo3 = pd.concat([combo2, values], axis=1)

        combo3.to_file(land_use_mataura_gpkg)

        # # Remove geometry for now, but I might need it later...
        # combo4 = pd.DataFrame(combo3.drop('geometry', axis=1))
        # combo5 = pd.concat([combo4, awm])
        # combo5 = combo5.loc[~combo5.index.duplicated()]

        # combo_dict[catch_id] = combo5


























































