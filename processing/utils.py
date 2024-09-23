#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 09:19:21 2022

@author: mike
"""
import warnings
warnings.simplefilter(action='ignore', category=RuntimeWarning)

import copy
import pathlib
import os
import pickle
import io
from shapely.ops import unary_union
from shapely import intersection, difference, intersects, symmetric_difference, make_valid
import geobuf
import base64
import orjson
import zstandard as zstd
import geopandas as gpd
import pandas as pd
import numpy as np
import xarray as xr
import hdf5tools
from scipy import stats
import booklet
import tempfile
import geopandas as gpd

import params

###############################################
### Lake equations


def est_b(t):
    """

    """
    return 1 + 0.44*(t**0.13)


# def est_log_c_lake_p(c_p_in, t, z_max):
#     """

#     """
#     if z_max > 7.5:
#         return np.log10(c_p_in)/est_b(t)
#     else:
#         return np.log10(c_p_in)


# def est_log_c_lake_n(c_n_in, z_max):
#     """

#     """
#     return 1.6 + 0.54*np.log10(c_n_in) - 0.41*np.log10(z_max)


# def est_log_c_lake_chla(log_c_lake_p, log_c_lake_n):
#     """

#     """
#     return -1.8 + 0.7*log_c_lake_n + 0.55*log_c_lake_p


# def est_d_lake_secchi(log_c_lake_chla, z_max, u, fetch):
#     """

#     """
#     if z_max >= 20:
#         return (3.46 - 1.53*log_c_lake_chla)**2
#     else:
#         return (3.46 - 0.74*log_c_lake_chla - 0.35*np.log10((fetch*(u**2))/z_max))**2


def est_r_lake_p(r_p_in, t, z_max):
    """

    """
    if z_max > 7.5:
        return r_p_in**(1/est_b(t))
    else:
        return r_p_in


def est_r_lake_n(r_n_in):
    """

    """
    return r_n_in**0.54


def est_r_lake_chla(r_p_lake, r_n_lake):
    """

    """
    return (r_p_lake**0.55) * (r_n_lake**0.7)


def est_d_lake_secchi_scenario(r_chla_lake, d_secchi_lake_current, z_max):
    """

    """
    if z_max >= 20:
        return (np.log10(r_chla_lake**-1.53) + (d_secchi_lake_current**0.5))**2
    else:
        return (np.log10(r_chla_lake**-0.74) + (d_secchi_lake_current**0.5))**2


### Regional equations

def est_b_canterbury(t, z_max):
    """

    """
    if z_max > 7.5:
        b = 1 + 0.91888*(t**0.0205)
    else:
        b = 1 + 0.09288*(t**0.0205)

    return b


# def est_log_c_lake_p_waikato(c_p_in, t):
#     """

#     """
#     b = 1 + t**0.5

#     return 0.9217 + 0.6172 * (np.log10(c_p_in)/b)


# def est_log_c_lake_n_waikato(c_n_in, t):
#     """

#     """
#     b = 1 + t**0.5

#     return 2.3969 + 0.3564 * (np.log10(c_n_in)/b)


# def est_log_c_lake_p_canterbury(c_p_in, t, z_max):
#     """

#     """
#     b = est_b_canterbury(t, z_max)

#     return np.log10(c_p_in)/b


def est_r_lake_p_waikato(r_p_in, t):
    """

    """
    b = 1 + t**0.5

    return r_p_in**(0.6172/b)


def est_r_lake_n_waikato(r_n_in, t):
    """

    """
    b = 1 + t**0.5

    return r_n_in**(0.3564/b)


def est_r_lake_p_canterbury(r_p_in, t, z_max):
    """

    """
    b = est_b_canterbury(t, z_max)

    return r_p_in**(1/b)


def est_ind_scenario_conc(current_concs, conc_factors, region, t, z_max):
    """

    """
    ## Turn off regional estimates...because they make no sense...
    # region = ''

    ## Nitrogen
    r_n_in = conc_factors['nitrogen']
    c_n_lake_current = current_concs['TN']
    if region == 'Waikato':
        r_n_lake = est_r_lake_n_waikato(r_n_in, t)
    else:
        r_n_lake = est_r_lake_n(r_n_in)

    c_n_lake_scenario = c_n_lake_current * r_n_lake

    ## TP
    r_p_in = conc_factors['phosphorus']
    c_p_lake_current = current_concs['TP']
    if region == 'Waikato':
        r_p_lake = est_r_lake_p_waikato(r_p_in, t)
    elif region == 'Canterbury':
        r_p_lake = est_r_lake_p_canterbury(r_p_in, t, z_max)
    else:
        r_p_lake = est_r_lake_p(r_p_in, t, z_max)

    c_p_lake_scenario = c_p_lake_current * r_p_lake

    ## Chla
    c_chla_lake_current = current_concs['CHLA']

    r_chla_lake = est_r_lake_chla(r_p_lake, r_n_lake)

    c_chla_lake_scenario = c_chla_lake_current * r_chla_lake

    ## Secchi
    d_secchi_lake_current = current_concs['Secchi']

    d_secchi_lake_scenario = est_d_lake_secchi_scenario(r_chla_lake, d_secchi_lake_current, z_max)

    ## Package up the results
    res = {'TN': c_n_lake_scenario, 'TP': c_p_lake_scenario, 'CHLA': c_chla_lake_scenario, 'Secchi': d_secchi_lake_scenario}

    return res

#############################################
### Functions


def find_downstream(way_id, node_way, way):
    """

    """
    way_down = way_id
    ways_down = [way_down]
    append = ways_down.append

    old_down_node = -1

    while True:
        down_node = way[int(way_down)][-1]
        if old_down_node == down_node:
            break

        down_ways = list(node_way[int(down_node)])
        down_ways.remove(way_down)

        if down_ways:
            for down_way in down_ways:
                up_node = way[int(down_way)][0]
                if up_node == down_node:
                    append(down_way)
                    way_down = down_way
                    break
                else:
                    old_down_node = down_node
        else:
            break

    return ways_down


def area_weighted_mean(data, cat_col, cols):
    """

    """
    if isinstance(data, gpd.GeoDataFrame):
        if data.geom_type.isin(['MultiPolygon', 'Polygon']).all():
            tot_areas = data.geometry.area.groupby(data[cat_col]).transform('sum')
            ratio_areas = data.geometry.area/tot_areas
            new_data = data[cols].mul(ratio_areas.values, axis=0)
            new_data_agg = new_data.groupby(data[cat_col]).sum()

            return new_data_agg

        else:
            raise ValueError('Not all geometries are Polygons!')
    else:
        raise ValueError('data needs to be a geodataframe!')


def calc_stats(data):
    """

    """
    min1 = np.min(data)
    median1 = np.median(data)
    q25 = np.percentile(data, 25)
    q75 = np.percentile(data, 75)
    q5 = np.percentile(data, 5)
    q95 = np.percentile(data, 95)
    max1 = np.max(data)
    count1 = len(data)
    from_date = data.index.min().isoformat(sep=' ', timespec='minutes')
    to_date = data.index.max().isoformat(sep=' ', timespec='minutes')

    return {'min': min1, 'Q5': q5, 'Q25': q25, 'median': median1, 'Q75': q75, 'Q95': q95, 'max': max1, 'count': count1, 'from_date': from_date, 'to_date': to_date}


def gpd_to_feather(gdf, output):
    """

    """
    gdf.to_feather(output, compression='zstd', compression_level=1)


def df_to_feather(df, output):
    """

    """
    df.to_feather(output, compression='zstd', compression_level=1)


def read_pkl_zstd(obj, unpickle=False):
    """
    Deserializer from a pickled object compressed with zstandard.

    Parameters
    ----------
    obj : bytes or str
        Either a bytes object that has been pickled and compressed or a str path to the file object.
    unpickle : bool
        Should the bytes object be unpickled or left as bytes?

    Returns
    -------
    Python object
    """
    if isinstance(obj, (str, pathlib.Path)):
        with open(obj, 'rb') as p:
            dctx = zstd.ZstdDecompressor()
            with dctx.stream_reader(p) as reader:
                obj1 = reader.read()

    elif isinstance(obj, bytes):
        dctx = zstd.ZstdDecompressor()
        obj1 = dctx.decompress(obj)
    else:
        raise TypeError('obj must either be a str path or a bytes object')

    if unpickle:
        obj1 = pickle.loads(obj1)

    return obj1


def write_pkl_zstd(obj, file_path=None, compress_level=1, pkl_protocol=pickle.HIGHEST_PROTOCOL, retries=3):
    """
    Serializer using pickle and zstandard. Converts any object that can be pickled to a binary object then compresses it using zstandard. Optionally saves the object to disk. If obj is bytes, then it will only be compressed without pickling.

    Parameters
    ----------
    obj : any
        Any pickleable object.
    file_path : None or str
        Either None to return the bytes object or a str path to save it to disk.
    compress_level : int
        zstandard compression level.

    Returns
    -------
    If file_path is None, then it returns the byte object, else None.
    """
    if isinstance(obj, bytes):
        p_obj = obj
    else:
        p_obj = pickle.dumps(obj, protocol=pkl_protocol)

    if isinstance(file_path, (str, pathlib.Path)):

        with open(file_path, 'wb') as f:
            cctx = zstd.ZstdCompressor(level=compress_level, write_content_size=True)
            with cctx.stream_writer(f, size=len(p_obj)) as compressor:
                compressor.write(p_obj)
    else:
        cctx = zstd.ZstdCompressor(level=compress_level, write_content_size=True)
        c_obj = cctx.compress(p_obj)

        return c_obj


def xr_concat(datasets):
    """
    A much more efficient concat/combine of xarray datasets. It's also much safer on memory.
    """
    # Get variables for the creation of blank dataset
    coords_list = []
    chunk_dict = {}

    for chunk in datasets:
        coords_list.append(chunk.coords.to_dataset())
        for var in chunk.data_vars:
            if var not in chunk_dict:
                dims = tuple(chunk[var].dims)
                enc = chunk[var].encoding.copy()
                dtype = chunk[var].dtype
                _ = [enc.pop(d) for d in ['original_shape', 'source'] if d in enc]
                var_dict = {'dims': dims, 'enc': enc, 'dtype': dtype, 'attrs': chunk[var].attrs}
                chunk_dict[var] = var_dict

    try:
        xr3 = xr.combine_by_coords(coords_list, compat='override', data_vars='minimal', coords='all', combine_attrs='override')
    except:
        xr3 = xr.merge(coords_list, compat='override', combine_attrs='override')

    # Create the blank dataset
    for var, var_dict in chunk_dict.items():
        dims = var_dict['dims']
        shape = tuple(xr3[c].shape[0] for c in dims)
        xr3[var] = (dims, np.full(shape, np.nan, var_dict['dtype']))
        xr3[var].attrs = var_dict['attrs']
        xr3[var].encoding = var_dict['enc']

    # Update the attributes in the coords from the first ds
    for coord in xr3.coords:
        xr3[coord].encoding = datasets[0][coord].encoding
        xr3[coord].attrs = datasets[0][coord].attrs

    # Fill the dataset with data
    for chunk in datasets:
        for var in chunk.data_vars:
            if isinstance(chunk[var].variable._data, np.ndarray):
                xr3[var].loc[chunk[var].transpose(*chunk_dict[var]['dims']).coords.indexes] = chunk[var].transpose(*chunk_dict[var]['dims']).values
            elif isinstance(chunk[var].variable._data, xr.core.indexing.MemoryCachedArray):
                c1 = chunk[var].copy().load().transpose(*chunk_dict[var]['dims'])
                xr3[var].loc[c1.coords.indexes] = c1.values
                c1.close()
                del c1
            else:
                raise TypeError('Dataset data should be either an ndarray or a MemoryCachedArray.')

    return xr3


def get_directly_upstream_ways(way_id, node_way, way, way_index):
    """

    """
    ways_up = set([way_id])

    new_ways = set(way_index[int(way_id)]).difference(ways_up)

    down_node = way[int(way_id)][-1]
    down_ways = set(node_way[down_node])

    new_ways = new_ways.difference(down_ways)

    return new_ways


def discrete_resample(df, freq_code, agg_fun, remove_inter=False, **kwargs):
    """
    Function to properly set up a resampling class for discrete data. This assumes a linear interpolation between data points.

    Parameters
    ----------
    df: DataFrame or Series
        DataFrame or Series with a time index.
    freq_code: str
        Pandas frequency code. e.g. 'D'.
    agg_fun : str
        The aggregation function to be applied on the resampling object.
    **kwargs
        Any keyword args passed to Pandas resample.

    Returns
    -------
    Pandas DataFrame or Series
    """
    if isinstance(df, (pd.Series, pd.DataFrame)):
        if isinstance(df.index, pd.DatetimeIndex):
            reg1 = pd.date_range(df.index[0].ceil(freq_code), df.index[-1].floor(freq_code), freq=freq_code)
            reg2 = reg1[~reg1.isin(df.index)]
            if isinstance(df, pd.Series):
                s1 = pd.Series(np.nan, index=reg2)
            else:
                s1 = pd.DataFrame(np.nan, index=reg2, columns=df.columns)
            s2 = pd.concat([df, s1]).sort_index()
            s3 = s2.interpolate('time')
            s4 = (s3 + s3.shift(-1))/2
            s5 = s4.resample(freq_code, **kwargs).agg(agg_fun).dropna()

            if remove_inter:
                index1 = df.index.floor(freq_code).unique()
                s6 = s5[s5.index.isin(index1)].copy()
            else:
                s6 = s5
        else:
            raise ValueError('The index must be a datetimeindex')
    else:
        raise TypeError('The object must be either a DataFrame or a Series')

    return s6


def dtl_correction(data, dtl_method='trend'):
    """
    The method to use to convert values below a detection limit to numeric. Used for water quality results. Options are 'half' or 'trend'. 'half' simply halves the detection limit value, while 'trend' uses half the highest detection limit across the results when more than 40% of the values are below the detection limit. Otherwise it uses half the detection limit.
    """
    new_data_list = []
    append = new_data_list.append
    for i, df in data.groupby(['lawa_id', 'indicator']):
        if df.censor_code.isin(['greater_than', 'less_than']).any():
            greater1 = df.censor_code == 'greater_than'
            df.loc[greater1, 'value'] = df.loc[greater1, 'value'] * 1.5

            less1 = df.censor_code == 'less_than'
            if less1.sum() > 0:
                df.loc[less1, 'value'] = df.loc[less1, 'value'] * 0.5
                if dtl_method == 'trend':
                    df1 = df.loc[less1]
                    count1 = len(df)
                    count_dtl = len(df1)
                    # count_dtl_val = df1['value'].nunique()
                    dtl_ratio = np.round(count_dtl / float(count1), 2)
                    if dtl_ratio >= 0.4:
                        dtl_val = df1['value'].max()
                        df.loc[(df['value'] < dtl_val) | less1, 'value'] = dtl_val

        append(df)

    new_data = pd.concat(new_data_list)

    return new_data


def process_land_cover_per_catch(way_id, data):
    """

    """
    with booklet.open(params.sites_reach_mapping_path) as reaches:
        reaches0 = reaches[way_id]

    ## Break into 1st order catchments
    geos = data.geometry.values

    with booklet.open(params.sites_catch_minor_path) as catches_minor:
        minor_catch_list = []
        for reach in reaches0:
            poly1 = catches_minor[reach]
            gs = gpd.GeoSeries([poly1], crs=4326).to_crs(2193).iloc[0]
            intersect1 = intersection(geos, gs)
            areas = []
            for i in intersect1:
                if i is None:
                    areas.append(-1)
                else:
                    areas.append(i.area)
            # areas = np.asarray([i.area for i in intersect1])
            max_index = np.argmax(areas)
            catch_data = data.iloc[[max_index]].copy()
            catch_data['geometry'] = [gs]
            catch_data['nzsegment'] = [reach]
            minor_catch_list.append(catch_data)

    minor_catch0 = pd.concat(minor_catch_list)
    # minor_catch0['area_m2'] = minor_catch0.area.round().astype('int32')
    # minor_catch0['catch_id'] = way_id

    return minor_catch0


def calc_farm_yields(poly, snb_dairy1, farm_yields, slope_geo0, moisture_geo0):
    """

    """
    # SnB and dairy
    # snb_dairy1 = farm_geo.loc[farm_geo.sindex.query(poly, predicate="intersects")].copy()
    # if not snb_dairy1.empty:
    #     slope_geo0 = slope_geo.loc[slope_geo.sindex.query(poly, predicate="intersects")]
    #     moisture_geo0 = moisture_geo.loc[moisture_geo.sindex.query(poly, predicate="intersects")]

    snb_dairy1b = intersection(snb_dairy1.geometry.tolist(), poly)
    snb_dairy1['geometry'] = snb_dairy1b
    snb_dairy1 = snb_dairy1[~snb_dairy1.geometry.is_empty].copy()
    snb_dairy2 = snb_dairy1.sjoin(slope_geo0, how='left', rsuffix='slope')
    snb_dairy3 = snb_dairy2.sjoin(moisture_geo0, how='left', rsuffix='moisture')
    snb_dairy4 = snb_dairy3.merge(farm_yields, on=['land_cover', 'slope', 'moisture']).drop(['index_slope', 'index_moisture', 'slope', 'moisture'], axis=1)
    snb_yields0 = area_weighted_mean(snb_dairy4, 'typology', ['phosphorus_yield', 'nitrogen_yield']).round(3)
    snb_dairy5 = snb_dairy4.dissolve('typology').drop(['phosphorus_yield', 'nitrogen_yield'], axis=1)
    snb_dairy6 = snb_dairy5.merge(snb_yields0, on='typology').reset_index()
    # snb_dairy1['geometry'] = snb_dairy1['geometry'].buffer(0.1).simplify(10).make_valid()
    snb_dairy6['geometry'] = snb_dairy6['geometry'].buffer(0).make_valid()

    return snb_dairy6


def calc_lcdb_yields(poly, lcdb_geo):
    """

    """
    # LCDB
    lcdb3 = lcdb_geo.loc[lcdb_geo.sindex.query(poly, predicate="intersects")].copy()
    if not lcdb3.empty:
        lcdb3b = intersection(lcdb3.geometry.tolist(), poly)
        lcdb3['geometry'] = lcdb3b
        lcdb3 = lcdb3[~lcdb3.geometry.is_empty].copy()
        lcdb3 = lcdb3.dissolve('typology').reset_index()
        # lcdb3['geometry'] = lcdb3.buffer(0.1).simplify(10).make_valid()
        lcdb3['geometry'] = lcdb3.buffer(0).make_valid()

    return lcdb3


def combine_yields(lake_id, poly):
    """

    """
    with booklet.open(params.lakes_lcdb_yields_blt_path) as f:
        lcdb3 = f[lake_id]
    with booklet.open(params.lakes_farm_yields_blt_path) as f:
        snb_dairy1 = f[lake_id]

    # Combo
    if (not snb_dairy1.empty) and (not lcdb3.empty):
        snb_dairy1['geometry'] = snb_dairy1.buffer(0.01).simplify(10).make_valid()
        lcdb3['geometry'] = lcdb3.buffer(0.01).simplify(10).make_valid()
        lcdb3 = lcdb3.overlay(snb_dairy1, how='difference', keep_geom_type=True)
        combo1 = pd.concat([snb_dairy1, lcdb3])
    elif snb_dairy1.empty:
        lcdb3['geometry'] = lcdb3.buffer(0.01).simplify(10).make_valid()
        combo1 = lcdb3
    else:
        snb_dairy1['geometry'] = snb_dairy1.buffer(0.01).simplify(10).make_valid()
        combo1 = snb_dairy1

    # Add in unimproveable areas
    diff_geo = difference(poly, make_valid(combo1.geometry.unary_union))
    if diff_geo is None:
        diff_gpd = gpd.GeoDataFrame([{'phosphorus_yield': 0.3, 'nitrogen_yield': 2}], geometry=[poly], crs=2193)
    else:
        diff_gpd = gpd.GeoDataFrame([{'phosphorus_yield': 0.3, 'nitrogen_yield': 2}], geometry=[diff_geo], crs=2193)
    diff_gpd['land_cover'] = ['Native Vegetation']
    diff_gpd['typology'] = ['Native Vegetation']
    combo5 = pd.concat([combo1.drop(['farm_type'], axis=1), diff_gpd])
    combo5['area_ha'] = (combo5.geometry.area*0.0001).round(2)

    # Update lcdb if possible
    lcs = combo5.land_cover.unique()
    if ('Sheep and Beef' in lcs) and ('Low Producing Grassland' in lcs):
        snb_mean = combo5.loc[combo5.land_cover == 'Sheep and Beef', ['nitrogen_yield', 'phosphorus_yield']].mean()
        combo5.loc[combo5.land_cover == 'Low Producing Grassland', ['nitrogen_yield', 'phosphorus_yield']] = snb_mean.values

    if ('Dairy' in lcs) and ('High Producing Exotic Grassland' in lcs):
        snb_mean = combo5.loc[combo5.land_cover == 'Dairy', ['nitrogen_yield', 'phosphorus_yield']].mean()
        combo5.loc[combo5.land_cover == 'High Producing Exotic Grassland', ['nitrogen_yield', 'phosphorus_yield']] = snb_mean.values

    # if catch_id in params.typo_corrections:
    #     typo_corr = params.typo_corrections[catch_id]
    #     for typo, corr in typo_corr.items():
    #         for ind, val in corr.items():
    #             combo5.loc[combo5.typology == typo, ind] = val

    return combo5


# def calc_yields(catch_id, poly, farm_geo, farm_yields, slope_geo, moisture_geo, lcdb_geo):
#     """

#     """
#     # SnB and dairy
#     snb_dairy1 = farm_geo.loc[farm_geo.sindex.query(poly, predicate="intersects")].copy()
#     if not snb_dairy1.empty:
#         snb_dairy1b = intersection(snb_dairy1.geometry.tolist(), poly)
#         snb_dairy1['geometry'] = snb_dairy1b
#         snb_dairy1 = snb_dairy1[~snb_dairy1.geometry.is_empty].copy()
#         snb_dairy2 = snb_dairy1.sjoin(slope_geo, how='left', rsuffix='slope')
#         snb_dairy3 = snb_dairy2.sjoin(moisture_geo, how='left', rsuffix='moisture')
#         snb_dairy4 = snb_dairy3.merge(farm_yields, on=['land_cover', 'slope', 'moisture']).drop(['index_slope', 'index_moisture', 'slope', 'moisture'], axis=1)
#         snb_yields0 = area_weighted_mean(snb_dairy4, 'typology', ['phosphorus_yield', 'nitrogen_yield']).round(3)
#         snb_dairy5 = snb_dairy4.dissolve('typology').drop(['phosphorus_yield', 'nitrogen_yield'], axis=1)
#         snb_dairy1 = snb_dairy5.merge(snb_yields0, on='typology').reset_index()
#         # snb_dairy1['geometry'] = snb_dairy1['geometry'].buffer(0.1).simplify(10).make_valid()
#         snb_dairy1['geometry'] = snb_dairy1['geometry'].buffer(0).make_valid()

#     # LCDB
#     lcdb3 = lcdb_geo.loc[lcdb_geo.sindex.query(poly, predicate="intersects")].copy()
#     if not lcdb3.empty:
#         lcdb3b = intersection(lcdb3.geometry.tolist(), poly)
#         lcdb3['geometry'] = lcdb3b
#         lcdb3 = lcdb3[~lcdb3.geometry.is_empty].copy()
#         lcdb3 = lcdb3.dissolve('typology').reset_index()
#         # lcdb3['geometry'] = lcdb3.buffer(0.1).simplify(10).make_valid()
#         lcdb3['geometry'] = lcdb3.buffer(0).make_valid()

#     # Combo
#     if (not snb_dairy1.empty) and (not lcdb3.empty):
#         try:
#             lcdb3 = lcdb3.overlay(snb_dairy1, how='difference', keep_geom_type=True)
#         except Exception as error:
#             print(catch_id)
#             raise error
#         combo1 = pd.concat([snb_dairy1, lcdb3])
#     elif snb_dairy1.empty:
#         combo1 = lcdb3
#     else:
#         combo1 = snb_dairy1

#     # Add in unimproveable areas
#     diff_geo = difference(poly, make_valid(combo1.geometry.unary_union))
#     if diff_geo is None:
#         diff_gpd = gpd.GeoDataFrame([{'phosphorus_yield': 0.3, 'nitrogen_yield': 2}], geometry=[poly], crs=2193)
#     else:
#         diff_gpd = gpd.GeoDataFrame([{'phosphorus_yield': 0.3, 'nitrogen_yield': 2}], geometry=[diff_geo], crs=2193)
#     diff_gpd['land_cover'] = ['Native Vegetation']
#     diff_gpd['typology'] = ['Native Vegetation']
#     combo5 = pd.concat([combo1.drop(['farm_type'], axis=1), diff_gpd])
#     combo5['area_ha'] = (combo5.geometry.area*0.0001).round(2)

#     # Update lcdb if possible
#     lcs = combo5.land_cover.unique()
#     if ('Sheep and Beef' in lcs) and ('Low Producing Grassland' in lcs):
#         snb_mean = combo5.loc[combo5.land_cover == 'Sheep and Beef', ['nitrogen_yield', 'phosphorus_yield']].mean()
#         combo5.loc[combo5.land_cover == 'Low Producing Grassland', ['nitrogen_yield', 'phosphorus_yield']] = snb_mean.values

#     if ('Dairy' in lcs) and ('High Producing Exotic Grassland' in lcs):
#         snb_mean = combo5.loc[combo5.land_cover == 'Dairy', ['nitrogen_yield', 'phosphorus_yield']].mean()
#         combo5.loc[combo5.land_cover == 'High Producing Exotic Grassland', ['nitrogen_yield', 'phosphorus_yield']] = snb_mean.values

#     # if catch_id in params.typo_corrections:
#     #     typo_corr = params.typo_corrections[catch_id]
#     #     for typo, corr in typo_corr.items():
#     #         for ind, val in corr.items():
#     #             combo5.loc[combo5.typology == typo, ind] = val

#     return combo5


def calc_reductions(lake_id, poly, lcdb1, snb_dairy1):
    """

    """
    # Land cover
    # lcdb1 = lcdb.loc[lcdb.sindex.query(poly, predicate="intersects")].copy()
    if not lcdb1.empty:
        lcdb1b = intersection(lcdb1.geometry.tolist(), poly)
        lcdb1['geometry'] = lcdb1b
        lcdb1 = lcdb1[~lcdb1.geometry.is_empty].copy()
        lcdb1 = lcdb1.dissolve('typology').reset_index()
        lcdb1['geometry'] = lcdb1.buffer(0.1).simplify(10).make_valid()

        # lc3['geometry'] = lc3.buffer(0.5).simplify(10)

    ## SnB and Dairy
    # snb_dairy1 = snb_dairy.loc[snb_dairy.sindex.query(poly, predicate="intersects")].copy()
    if not snb_dairy1.empty:
        snb_dairy1b = intersection(snb_dairy1.geometry.tolist(), poly)
        snb_dairy1['geometry'] = snb_dairy1b
        snb_dairy1 = snb_dairy1[~snb_dairy1.geometry.is_empty].copy()
        snb_dairy1 = snb_dairy1.dissolve('typology').reset_index()
        snb_dairy1['geometry'] = snb_dairy1['geometry'].buffer(0.1).simplify(10).make_valid()

    if (not snb_dairy1.empty) and (not lcdb1.empty):
        # diff_list = []
        # for geo1 in lcdb1.geometry:
        #     for geo2 in snb_dairy1.geometry:
        #         if intersects(geo1, geo2):
        #             geo3 = difference(geo1, geo2)
        #             diff_list.append(geo3)
        #         else:
        #             diff_list.append(geo1)
        lcdb1 = lcdb1.overlay(snb_dairy1, how='difference', keep_geom_type=True)
        combo2 = pd.concat([snb_dairy1, lcdb1])
    elif snb_dairy1.empty:
        combo2 = lcdb1
    else:
        combo2 = snb_dairy1

    # if catch_id in params.typo_red_corrections:
    #     typo_corr = params.typo_red_corrections[catch_id]
    #     for typo, corr in typo_corr.items():
    #         for ind, val in corr.items():
    #             combo2.loc[combo2.typology == typo, ind] = val

    ## Aggregate classes
    # combo3 = combo2.replace({'land_cover': params.combine_land_covers})
    combo4 = combo2.drop(['farm_type'], axis=1)
    diff_geo = difference(poly, make_valid(combo4.geometry.unary_union))
    if diff_geo is None:
        diff_gpd = gpd.GeoDataFrame([{'phosphorus_reduction': 0, 'nitrogen_reduction': 0}], geometry=[poly], crs=2193)
    else:
        diff_gpd = gpd.GeoDataFrame([{'phosphorus_reduction': 0, 'nitrogen_reduction': 0}], geometry=[diff_geo], crs=2193)
    diff_gpd['land_cover'] = ['Native Vegetation']
    diff_gpd['typology'] = ['Native Vegetation']
    combo5 = pd.concat([combo4, diff_gpd])
    combo5['area_ha'] = (combo5.geometry.area*0.0001).round(2)

    return combo5


def ecoli_spatial_yields(data, combo2, ecoli_factors):
    """

    """
    # with booklet.open(params.sites_yields_blt_path) as y:
    #     data = y[catch_id]
    data1 = data[['typology', 'geometry']].sjoin(combo2, how='left', predicate='intersects').drop(['index_right'], axis=1)
    missing = data1.loc[data1.elev.isnull(), ['typology', 'geometry']]
    if not missing.empty:
        miss1 = missing.sjoin_nearest(combo2, how='left', max_distance=300).drop(['index_right'], axis=1)
        data1 = pd.concat([data1, miss1])

    data2 = data1.groupby('typology')[['elev', 'drain']].median()
    data2['elev'] = pd.cut(data2['elev'], [-100, 350, 1000], labels=['LowElev', 'HighElev'])
    data2['drain'] = pd.cut(data2['drain'], [0, 3.5, 10], labels=['PoorlyDrained', 'WellDrained'])
    # data3 = pd.merge(data2.reset_index(), lu_map, on='lu_id').drop('lu_id', axis=1)

    data4 = data.drop(['phosphorus_yield', 'nitrogen_yield'], axis=1).merge(data2, on='typology').drop('typology', axis=1)

    data4['lu'] = data4.land_cover.replace(params.combine_land_covers)

    data5 = pd.merge(data4.drop('geometry', axis=1), ecoli_factors, on=['elev', 'drain', 'lu'], how='outer')
    data5.loc[data5.area_ha.isnull(), 'area_ha'] = 0
    # data5.loc[data5.area_ratio.isnull(), 'area_ratio'] = 0
    data5.loc[data5.ecoli_factor.isnull(), 'ecoli_factor'] = 0

    ## Add in simplified land use classes
    # data5['lu'] = data5.lu.replace(params.lu_to_land_cover)
    data5.loc[data5.land_cover.isnull(), 'land_cover'] = data5.loc[data5.land_cover.isnull(), 'lu']
    data5['land_use'] = data5.land_cover.replace(params.combine_land_covers)

    area_sum1 = data5.groupby(['land_cover', 'elev', 'drain'])['area_ha'].transform('sum')
    area_ratio1 = data5.area_ha/area_sum1
    area_ratio1.loc[area_ratio1.isnull()] = 1

    data5['ecoli_factor'] = area_ratio1 * data5['ecoli_factor']

    data6 = data5.groupby(['land_use', 'elev', 'drain'])[['area_ha', 'ecoli_factor']].sum()

    ## Add in reductions
    # data5['ecoli_reduction'] = data5.land_use.replace(params.ecoli_reductions)

    return data6





















































