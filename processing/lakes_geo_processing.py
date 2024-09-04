#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 16:17:00 2022

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
# import dbm
import booklet
# import shelve
import multiprocessing as mp
import concurrent.futures
import geobuf
from shapely.geometry import Point, Polygon, box, LineString, mapping, shape
import orjson

import utils, params

pd.options.display.max_columns = 10

##################################################
### preprocessing

cat_cols = ['Current5', 'GeomorphicType']
num_cols = ['MaxDepth', 'LakeArea', 'DecTemp', 'DecSolrad', 'Fetch', 'SumWind', 'CatBeech', 'CatGlacial', 'CatHard', 'CatPeat', 'CatPhos', 'CatSlope', 'CatAnnTemp', 'DirectDistCoast', 'ResidenceTime', 'Urban', 'Pasture', 'LakeElevation', 'MeanWind']
model_cols = cat_cols + num_cols


def lakes_points_poly_process():

    ## Lakes polygons
    lakes_poly0 = gpd.read_file(params.lakes_poly_path)[['LFENZID', 'Name', 'ResidenceTime', 'MaxDepth', 'LakeVolume', 'LakeAreaHa', 'geometry']].rename(columns={'Name': 'name', 'ResidenceTime': 'residence_time', 'MaxDepth': 'max_depth', 'LakeVolume': 'volume', 'LakeAreaHa': 'area_ha'}).copy()
    lakes_catch0 = gpd.read_file(params.lakes_catch_path)[['LID', 'catFlow']].rename(columns={'LID': 'LFENZID', 'catFlow': 'flow'})
    lakes_poly0 = lakes_poly0.merge(lakes_catch0, on='LFENZID')

    with booklet.open(params.lakes_reach_gbuf_path) as f:
        lake_ids = list(f.keys())

    lakes_poly0 = lakes_poly0[lakes_poly0.LFENZID.isin(lake_ids)].copy()
    lakes_poly0['LFENZID'] = lakes_poly0['LFENZID'].astype('int32')

    lakes_poly0['geometry'] = lakes_poly0.buffer(0).simplify(20).make_valid()
    lakes_poly0.loc[lakes_poly0.name == 'Lake Ototoa', 'LFENZID'] = 50270

    lakes_poly0 = lakes_poly0.dropna(subset=['LFENZID']).copy()
    lakes_poly0['residence_time'] = lakes_poly0['volume']/lakes_poly0['flow'] * 365
    # lakes_poly0.loc[lakes_poly0['residence_time'] < 1, 'residence_time'] = 1
    # lakes_poly0['residence_time'] = lakes_poly0['residence_time'].round().astype('int32')
    # lakes_poly0.loc[lakes_poly0['max_depth'].isnull(), 'max_depth'] = lakes_poly0['max_depth'].median()
    # lakes_poly0.loc[lakes_poly0['max_depth'] < 1, 'max_depth'] = 1
    # lakes_poly0['max_depth'] = lakes_poly0['max_depth'].round().astype('int32')
    lakes_poly0.loc[lakes_poly0.name.isnull(), 'name'] = 'No name'

    lakes_poly1 = lakes_poly0[lakes_poly0.residence_time.notnull() & lakes_poly0.max_depth.notnull()].copy()

    lakes_poly1['mean_depth'] = lakes_poly1.volume / lakes_poly1.area_ha / 10000
    lakes_poly1['p_residence_time'] = np.sqrt(lakes_poly1.residence_time)/(1 + np.sqrt(lakes_poly1.residence_time))
    lakes_poly1['n_residence_time'] = 1 - np.exp((-6.83 - lakes_poly1.residence_time)/(lakes_poly1.mean_depth))

    lakes_poly1 = lakes_poly1.drop_duplicates(subset=['LFENZID'])

    lakes_poly2 = lakes_poly1.loc[lakes_poly0.LFENZID.isin(lake_ids), ['LFENZID', 'name', 'residence_time', 'max_depth', 'mean_depth', 'p_residence_time', 'n_residence_time', 'geometry']].reset_index(drop=True).copy()

    with booklet.open(params.lakes_poly_gbuf_path, 'n', value_serializer='zstd', key_serializer='uint2', n_buckets=4001) as s:
        for LFENZID in lakes_poly2.LFENZID:
            geo = lakes_poly2[lakes_poly2.LFENZID == LFENZID].to_crs(4326).set_index('LFENZID', drop=False).__geo_interface__
            gbuf = geobuf.encode(geo)
            s[LFENZID] = gbuf

    lakes_poly2.to_file(params.lakes_poly_gpkg_path, index=False)

    ## Point locations of lakes
    # All lakes
    sites = lakes_poly2.copy()
    sites['geometry'] = sites.geometry.centroid.to_crs(4326)

    sites_geo = sites.set_index('LFENZID').__geo_interface__

    sites_gbuf = geobuf.encode(sites_geo)

    with open(params.lakes_points_gbuf_path, 'wb') as f:
        f.write(sites_gbuf)

    sites.to_file(params.lakes_points_gpkg_path, index=False)

    ## Point locations of monitoring sites
    # stdev0 = pd.read_csv(utils.lakes_stdev_moni_path)
    # site_loc0 = pd.read_csv(utils.lakes_raw_moni_data_csv_path, usecols=['LawaSiteID', 'SiteID', 'LFENZID', 'Latitude', 'Longitude',]).rename(columns={'LawaSiteID': 'lawa_id', 'SiteID': 'site_id', 'Latitude': 'lat', 'Longitude': 'lon'})
    # site_loc1 = site_loc0.drop_duplicates(subset=['lawa_id'])
    # site_loc2 = site_loc1[site_loc1.lawa_id.isin(stdev0.lawa_id.unique())].copy()
    # site_loc2['LFENZID'] = site_loc2['LFENZID'].astype('int32')

    # site_loc3 = vector.xy_to_gpd(['lawa_id', 'site_id', 'LFENZID'], 'lon', 'lat', site_loc2, 4326)

    # # All lakes
    # with booklet.open(utils.lakes_moni_sites_gbuf_path, 'n', key_serializer='uint4', value_serializer='zstd', n_buckets=4001) as f:
    #     for LFENZID in lakes_poly1.LFENZID:
    #         lake_geo = lakes_poly1[lakes_poly1.LFENZID == LFENZID].to_crs(4326).iloc[0].geometry
    #         site_loc4 = site_loc3[site_loc3.within(lake_geo)].set_index('site_id', drop=False).rename(columns={'site_id': 'tooltip'}).__geo_interface__
    #         gbuf = geobuf.encode(site_loc4)
    #         f[LFENZID] = gbuf

    # site_loc3.to_file(utils.lakes_moni_sites_gpkg_path)





















































