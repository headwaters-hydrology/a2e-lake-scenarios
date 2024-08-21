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



def delin_lakes_catch_rec():
    # break_points = gpd.read_file(utils.catch_break_points_gpkg).to_crs(4326)
    w0 = nzrec.Water(params.nzrec_data_path)
    rec_rivers0 = gpd.read_feather(params.rec_rivers_feather)

    stream_orders = {way_id: v['Strahler stream order'] for way_id, v in w0._way_tag.items()}

    ways_3rd_up = set([i for i, v in stream_orders.items() if v > 2])

    way = {k: v for k, v in w0._way.items()}
    way_index = {k: v for k, v in w0._way_index.items()}
    node_way = {k: v for k, v in w0._node_way_index.items()}
    # catch = {k: v for k, v in w0._catch.items()}
    catch = w0._catch

    # wq_sites = gpd.read_file(params.wq_sites_gpkg_path)

    # end_segs = wq_sites.nzsegment.unique()

    lakes_poly = gpd.read_file(params.lakes_poly_path)

    ## Lake inflows and outflows
    inflows = {}
    outflows = {}
    all_segs = {}
    for LFENZID in lakes_poly.LFENZID:
        if LFENZID > 0:
            LFENZID = int(LFENZID)
            print(LFENZID)

            lake_poly_geo = lakes_poly[lakes_poly.LFENZID == LFENZID]
            geom = lake_poly_geo.iloc[0]['geometry']

            # Inflows
            query_geom_ids = set(rec_rivers0[rec_rivers0.intersects(geom.boundary)].nzsegment.tolist())

            if query_geom_ids:
                inflows0 = set()
                if len(query_geom_ids) == 1:
                    inflows0.update(query_geom_ids)
                else:
                    geom_lat_lon = lake_poly_geo.to_crs(4326).iloc[0]['geometry']
                    for way_id in query_geom_ids:
                        nodes = way[way_id]
                        up_node = nodes[0]
                        up_node_coords = Point(w0._node[up_node]*0.0000001)
                        down_node = nodes[-1]
                        down_node_coords = Point(w0._node[down_node]*0.0000001)
    
                        if not shapely.intersects(geom_lat_lon, down_node_coords):
                            continue
                        if shapely.intersects(geom_lat_lon, up_node_coords):
                            continue
    
                        inflows0.add(way_id)

                    # Exceptions for small lakes
                    if geom.area < 500000:
                        diff_ways = query_geom_ids.difference(inflows0)
                        # diff_ways_inflows = set()
                        for way_id in diff_ways:
                            down_ways = utils.find_downstream(way_id, node_way, way)
                            if not inflows0.intersection(down_ways):
                                nodes = way[way_id]
                                up_node = nodes[0]
                                up_node_coords = Point(w0._node[up_node]*0.0000001)
                                down_node = nodes[-1]
                                down_node_coords = Point(w0._node[down_node]*0.0000001)
            
                                if not shapely.intersects(geom_lat_lon, down_node_coords) and not shapely.intersects(geom_lat_lon, up_node_coords):
                                    inflows0.add(way_id)

                    # Check for segments upstream of other inflows
                    for way_id in copy(inflows0):
                        up_ways = nzrec.utils.find_upstream(way_id, node_way, way, way_index)
                        up_ways.remove(way_id)
                        combo = inflows0.intersection(up_ways)
                        if combo:
                            for i in combo:
                                inflows0.remove(i)
                        # intersect0 = query_geom_ids.intersection(up_ways)
                        # if len(intersect0) == 1:
                        #     inflows0.add(way_id)

                    # if not inflows0:
                    #     if LFENZID != 15034:
                    #         d
                        # up_ways = nzrec.utils.find_upstream(way_id, node_way, way, way_index)

                inflows[LFENZID] = inflows0

            # Outflows
            query_geom_ids = set(rec_rivers0[rec_rivers0.intersects(geom)].nzsegment.tolist())

            outflows0 = set()
            if query_geom_ids:
                if len(query_geom_ids) == 1:
                    outflows0.update(query_geom_ids)
                else:
                    for way_id in query_geom_ids:
                        down_ways = utils.find_downstream(way_id, node_way, way)
                        if len(down_ways) > 1:
                            way_id_outflows = set()
                            for down_way in reversed(down_ways[1:]):
                                if down_way in query_geom_ids:
                                    way_id_outflows.add(down_way)
                                    break
                            if not way_id_outflows:
                                outflows0.add(way_id)
                            else:
                                outflows0.update(way_id_outflows)
                        else:
                            outflows0.add(way_id)

                outflows[LFENZID] = outflows0

                if not inflows[LFENZID]:
                    inflows[LFENZID].update(outflows0)

            # All segs in catchment
            if outflows0:
                all_segs0 = set()
                for way_id in outflows0:
                    up_ways = nzrec.utils.find_upstream(way_id, node_way, way, way_index)
                    all_segs0.update(up_ways)

                all_segs[LFENZID] = all_segs0


    ## Save results
    rec_rivers1 = rec_rivers0.set_index('nzsegment')

    geo_inflows_list  = []
    geo_outflows_list  = []
    # geo_all_list  = []
    for LFENZID in inflows:
        # Inflows
        inflows0 = inflows[LFENZID]
        rec0 = rec_rivers1.loc[list(inflows0)].copy()
        rec0['LFENZID'] = LFENZID
        geo_inflows_list.append(rec0)

        # Outflows
        outflows0 = outflows[LFENZID]
        rec0 = rec_rivers1.loc[list(outflows0)].copy()
        rec0['LFENZID'] = LFENZID
        geo_outflows_list.append(rec0)

        # All
        # all_segs0 = all_segs[LFENZID]
        # rec0 = rec_rivers1.loc[list(all_segs0)].copy()
        # rec0['LFENZID'] = LFENZID
        # geo_all_list.append(rec0)

    geo_inflows = pd.concat(geo_inflows_list)
    geo_outflows = pd.concat(geo_outflows_list)
    # geo_all = pd.concat(geo_all_list)

    geo_inflows.to_file(params.lake_inflows_gpkg)
    geo_outflows.to_file(params.lake_outflows_gpkg)
    # geo_all.to_file(params.lake_catch_segs_shp)

    # Booklet
    with booklet.open(params.lake_inflows_blt, 'n', value_serializer='pickle_zstd', key_serializer='uint4', n_buckets=1607) as reaches:
        for LFENZID, branches in inflows.items():
            reaches[LFENZID] = branches

    with booklet.open(params.lake_outflows_blt, 'n', value_serializer='pickle_zstd', key_serializer='uint4', n_buckets=1607) as reaches:
        for LFENZID, branches in outflows.items():
            reaches[LFENZID] = branches

    with booklet.open(params.lakes_reach_mapping_path, 'n', value_serializer='pickle_zstd', key_serializer='uint4', n_buckets=1607) as reaches:
        for LFENZID, branches in all_segs.items():
            reaches[LFENZID] = branches


    ## Catchment Delineation
    reach_gbuf_dict = {}
    catches_major_dict = {}
    # catches_major_4th_dict = {}
    catches_minor_dict = {}
    for LFENZID, ways_up in all_segs.items():
        print(LFENZID)

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

        reach_gbuf_dict[LFENZID] = gbuf

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
        catches_major_dict[LFENZID] = geo


    # Reach geobufs in blt
    with booklet.open(params.lakes_reach_gbuf_path, 'n', key_serializer='uint4', value_serializer='zstd', n_buckets=1607) as f:
        for LFENZID, gbuf in reach_gbuf_dict.items():
            f[LFENZID] = gbuf

    # Catchment gpds in blt
    with booklet.open(params.lakes_catch_major_path, 'n', key_serializer='uint4', value_serializer='wkb_zstd', n_buckets=1607) as f:
        for LFENZID, catches in catches_major_dict.items():
            f[LFENZID] = catches

    # Catchments geobuf
    with booklet.open(params.lakes_catch_major_gbuf_path, 'n', key_serializer='uint4', value_serializer='zstd', n_buckets=1607) as f:
        for LFENZID, catches in catches_major_dict.items():
            gdf = gpd.GeoDataFrame([{'LFENZID': way_id}], geometry=[catches], crs=4326).set_index('LFENZID', drop=False)
            gjson = orjson.loads(gdf.to_json())
            gbuf = geobuf.encode(gjson)
            f[LFENZID] = gbuf

    catch_ids = list(catches_major_dict.keys())
    rec_shed = gpd.GeoDataFrame(catch_ids, geometry=list(catches_major_dict.values()), crs=4326, columns=['LFENZID'])
    # rec_shed['geometry'] = rec_shed.simplify(0.0004)

    rec_shed.to_file(params.lakes_catch_gpkg_path)

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
    with booklet.open(params.lakes_catch_minor_path, 'n', key_serializer='uint4', value_serializer='pickle_zstd', n_buckets=500007) as f:
        for LFENZID, branches in catches_minor_dict.items():
            f[LFENZID] = branches


































































