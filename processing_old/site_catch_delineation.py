#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 08:19:06 2024

@author: mike
"""
import os
from gistools import vector, rec
import geopandas as gpd
import pandas as pd
import numpy as np
from copy import copy
import nzrec
import booklet
import orjson
import geobuf
import shapely
from shapely.geometry import Point, Polygon, box, LineString, mapping

import utils, params

pd.options.display.max_columns = 10


#############################################
### Rivers



def delin_sites_catch_rec():
    # break_points = gpd.read_file(utils.catch_break_points_gpkg).to_crs(4326)
    w0 = nzrec.Water(params.nzrec_data_path)

    stream_orders = {way_id: v['Strahler stream order'] for way_id, v in w0._way_tag.items()}

    ways_3rd_up = set([i for i, v in stream_orders.items() if v > 2])

    way = {k: v for k, v in w0._way.items()}
    # way_index = {k: v for k, v in w0._way_index.items()}
    # node_way = {k: v for k, v in w0._node_way_index.items()}
    catch = {k: v for k, v in w0._catch.items()}

    wq_sites = gpd.read_file(params.wq_sites_gpkg_path)

    end_segs = wq_sites.nzsegment.unique()

    ## Get the segs associated with the break points and remove he end ones from the end_segs
    # break_segs = set()
    # for i, row in break_points.iterrows():
    #     coords = np.array(row.geometry.coords)[0]
    #     way1 = w0.nearest_way(coords)
    #     break_segs.update(way1.ways)

    # for seg in break_segs:
    #     if seg in end_segs:
    #         end_segs.remove(seg)

    # ways_1st_2nd = set([i for i, v in stream_orders.items() if v < 3])

    # way_index_3rd_up = {way_id: set(v).intersection(ways_3rd_up) for way_id, v in way_index.items() if way_id in ways_3rd_up}

    ## Associate the 1st and 2nd order streams with the downstream 3rd order or greater streams
    # up_3rd_reaches_dict = {}
    # for way_id in ways_1st_2nd:
    #     new_ways = way_index[way_id]

    #     for w in new_ways:
    #         so = stream_orders[w]
    #         if so > 2:
    #             ways_up = nzrec.utils.find_upstream(way_id, node_way, way, way_index)
    #             if w in up_3rd_reaches_dict:
    #                 up_3rd_reaches_dict[w].update(ways_up)
    #             else:
    #                 up_3rd_reaches_dict[w] = ways_up
    #             break

    ## Aggregate the catchments and reassign
    # catch_aggs = {}
    # for way_id, branches in up_3rd_reaches_dict.items():
    #     branches.add(way_id)

    #     if way_id in catch_aggs:
    #         catch_aggs[way_id].update(branches)
    #     else:
    #         catch_aggs[way_id] = branches

    # for way_id, branches in catch_aggs.items():
    #     geo = shapely.ops.unary_union([catch[i] for i in branches])
    #     catch[way_id] = geo

    ## Delineate and aggregate 1st and 2nd order streams to the 3rd
    # reaches_dict = {}
    # missing_set = set()
    # for way_id in end_segs:
    #     # print(way_id)

    #     if way_id in ways_3rd_up:

    #         way0 = w0.add_way(way_id)
    #         all_up_ways = way0.upstream().ways

    #         reaches_dict[int(way_id)] = all_up_ways
    #     else:
    #         missing_set.add(way_id)

    reaches_dict = {}
    for way_id in end_segs:
        # print(way_id)

        way0 = w0.add_way(way_id)
        all_up_ways = way0.upstream().ways

        reaches_dict[int(way_id)] = all_up_ways


    with booklet.open(params.sites_reach_mapping_path, 'n', value_serializer='pickle_zstd', key_serializer='uint4', n_buckets=1607) as reaches:
        for way_id, branches in reaches_dict.items():
            reaches[way_id] = branches

    ## Delineate overall catchments
    reach_gbuf_dict = {}
    catches_major_dict = {}
    # catches_major_4th_dict = {}
    catches_minor_dict = {}
    for way_id, ways_up in reaches_dict.items():
        # print(way_id)

        # Reaches - only 3rd order and greater for rendering
        geo = []
        data = []
        for w in ways_up:
            if w in ways_3rd_up:
                nodes = way[w]
                geo.append(LineString(np.array([w0._node[i] * 0.0000001 for i in nodes])).simplify(0.0004))
                data.append({'nzsegment': int(w)})
        # data = [{'nzsegment': int(i)} for i in ways_up]

        if data:
            gdf = gpd.GeoDataFrame(data, geometry=geo, crs=4326).set_index('nzsegment', drop=False)
            gjson = orjson.loads(gdf.to_json())
            gbuf = geobuf.encode(gjson)
        else:
            gbuf = b''

        reach_gbuf_dict[way_id] = gbuf

        # Catchments
        geos = []
        append = geos.append
        for i in ways_up:
            geo0 = catch[i]
            if not geo0.is_valid:
                geo0 = shapely.make_valid(geo0)
            append(geo0)
            if i not in catches_minor_dict:
                catches_minor_dict[i] = geo0


        # geos = [catch[i] for i in ways_up]
        # catch1 = gpd.GeoDataFrame(ways_up, geometry=geos, crs=4326, columns=['nzsegment'])
        # catches_minor_dict[way_id] = catch1

        geo = shapely.ops.unary_union(geos).buffer(0.00001).simplify(0.0004)
        catches_major_dict[way_id] = geo

        # if stream_orders[way_id] > 3:
        #     catches_major_4th_dict[way_id] = geo


    # Reach geobufs in blt
    with booklet.open(params.sites_reach_gbuf_path, 'n', key_serializer='uint4', value_serializer='zstd', n_buckets=1607) as f:
        for way_id, gbuf in reach_gbuf_dict.items():
            f[way_id] = gbuf

    # Catchment gpds in blt
    with booklet.open(params.sites_catch_major_path, 'n', key_serializer='uint4', value_serializer='wkb_zstd', n_buckets=1607) as f:
        for way_id, catches in catches_major_dict.items():
            f[way_id] = catches

    # Catchments geobuf
    with booklet.open(params.sites_catch_major_gbuf_path, 'n', key_serializer='uint4', value_serializer='zstd', n_buckets=1607) as f:
        for way_id, catches in catches_major_dict.items():
            gdf = gpd.GeoDataFrame([{'nzsegment': way_id}], geometry=[catches], crs=4326).set_index('nzsegment', drop=False)
            gjson = orjson.loads(gdf.to_json())
            gbuf = geobuf.encode(gjson)
            f[way_id] = gbuf

    catch_ids = list(catches_major_dict.keys())
    rec_shed = gpd.GeoDataFrame(catch_ids, geometry=list(catches_major_dict.values()), crs=4326, columns=['nzsegment'])
    # rec_shed['geometry'] = rec_shed.simplify(0.0004)

    rec_shed.to_file(params.sites_catch_gpkg_path)

    # gjson = orjson.loads(rec_shed.set_index('nzsegment').to_json())

    # with open(params.sites_catch_major_gbuf_path, 'wb') as f:
    #     f.write(geobuf.encode(gjson))

    # catch_ids = list(catches_major_4th_dict.keys())
    # rec_shed = gpd.GeoDataFrame(catch_ids, geometry=list(catches_major_4th_dict.values()), crs=4326, columns=['nzsegment'])
    # rec_shed['geometry'] = rec_shed.simplify(0.0004)

    # gjson = orjson.loads(rec_shed.set_index('nzsegment').to_json())

    # with open(utils.assets_path.joinpath('rivers_catchments_4th.pbf'), 'wb') as f:
    #     f.write(geobuf.encode(gjson))

    ## Produce a file grouped by all catchments as geodataframes
    with booklet.open(params.sites_catch_minor_path, 'n', key_serializer='uint4', value_serializer='pickle_zstd', n_buckets=500007) as f:
        for way_id, branches in catches_minor_dict.items():
            f[way_id] = branches


































































