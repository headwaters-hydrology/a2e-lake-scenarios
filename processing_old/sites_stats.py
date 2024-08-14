#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 09:51:53 2024

@author: mike
"""
import os
import geopandas as gpd
import pandas as pd
import numpy as np
import booklet
import npsfm

import utils, params

pd.options.display.max_columns = 10


#############################################
### Rivers


def wq_data_stats():
    """

    """
    # nps = npsfm.NPSFM(params.data_path, download_files=True)

    # wq_sites = gpd.read_file(params.wq_sites_gpkg_path)
    # wq_sites = wq_sites.set_index('LawaSiteID')

    sites_stats = {}
    with booklet.open(params.wq_data_app_path, 'n', key_serializer='orjson', value_serializer='pickle_zstd', n_buckets=4391) as write:
        with booklet.open(params.wq_data_path) as f:
            for site_id, data in f.items():
                # nzsegment = wq_sites.loc[site_id, 'nzsegment']
                for grp, ts_data in data.groupby(['LawaSiteID', 'Indicator']):
                    ts_data0 = ts_data.set_index('SampleDateTime')['value']
                    stats = utils.calc_stats(ts_data0.dropna())
    
                    # if grp[1] in params.nps_maping:
                    #     parameter = params.nps_maping[grp[1]]
                    #     _ = nps.add_limits('river', parameter, nzsegment)
                    #     stats2 = nps.add_stats(ts_data0)
                    #     stats.update(stats2)
    
                    sites_stats[grp] = stats
                    write[grp] = ts_data0


    with booklet.open(params.wq_data_stats_path, 'n', key_serializer='orjson', value_serializer='orjson', n_buckets=4391) as f:
        for grp, stats in sites_stats.items():
            f[grp] = stats













































