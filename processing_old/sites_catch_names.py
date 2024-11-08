#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 08:51:20 2022

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

# name_exceptions = {3076139: '3076139 - Pokaiwhenua Stream'}


def sites_assign_catch_names():
    """

    """
    w0 = nzrec.Water(params.nzrec_data_path)

    with booklet.open(params.sites_catch_major_path) as f:
        catch_ids = set(f.keys())

    catch_names = {}
    for catch_id in catch_ids:
        # catch_name = str(catch_id)
        catch_name = 'Unnamed'
        tags = w0._way_tag[catch_id]
        if 'Catchment name' in tags:
            name = tags['Catchment name']
            if name is not None:
                if 'Unnamed' not in name:
                    catch_name = name
        catch_names[catch_id] = catch_name

    # for catch_id, name in name_exceptions.items():
    #     catch_names[catch_id] = name

    with booklet.open(params.sites_catch_names_path, 'n', key_serializer='uint4', value_serializer='str', n_buckets=1607) as f:
        for catch_id, catch_name in catch_names.items():
            f[catch_id] = catch_name
















































