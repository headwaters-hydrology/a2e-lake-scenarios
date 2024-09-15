#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 17:23:53 2023

@author: mike
"""
import os
import xarray as xr
from gistools import vector, rec
import geopandas as gpd
import pandas as pd
import numpy as np
import hdf5tools
import booklet
import nzrec

import utils, params

pd.options.display.max_columns = 10

#################################################
### Process loads


# segs = [3087421, 3087553, 3088075, 3088017, 3088076, 3095138]
# way_id = 11018765

tags = ['Climate class', 'Geology class', 'Topography class']


def process_reference_conc_rec():
    """

    """
    ## Determine all the segs
    # catches_minor = booklet.open(params.sites_catch_minor_path)
    # segs = np.asarray(list(catches_minor.keys()))
    # catches_minor.close()

    ## Import flows for load calcs
    # flows_list = []
    # append = flows_list.append
    # with booklet.open(params.river_flows_rec_diff_path) as f:
    #     for catch_id, flow in f.items():
    #         append((catch_id, flow))

    # flows = pd.DataFrame(flows_list, columns=['nzsegment', 'flow'])

    ### Reference conc - crazy...
    w0 = nzrec.Water(params.nzrec_data_path)

    # class_list = []
    # for seg in flows.nzsegment:
    #     data = w0._way_tag[seg]
    #     data1 = {'nzsegment': seg}
    #     data1.update({tag.split(' ')[0]: val for tag, val in data.items() if (tag in tags)})
    #     class_list.append(data1)

    class_list = []
    for seg, data in w0._way_tag.items():
        data1 = {'nzsegment': seg}
        data1.update({tag.split(' ')[0]: val for tag, val in data.items() if (tag in tags)})
        class_list.append(data1)

    seg_class0 = pd.DataFrame(class_list)
    seg_class0.loc[seg_class0['Climate'].isnull(), ['Climate', 'Topography', 'Geology']] = 'Other'

    ref_conc30 = pd.read_csv(params.rivers_ref_conc3_csv_path, header=[0, 1])

    ref_conc30['Climate'] = ref_conc30.REC.REC.str[:2]
    ref_conc30['Topography'] = ref_conc30.REC.REC.str[2:]

    ref_conc31 = ref_conc30.set_index([ref_conc30.param.param, ref_conc30['Climate'], ref_conc30['Topography']]).drop(['REC', 'param', 'Climate', 'Topography'], axis=1).stack()
    cols = list(ref_conc31.index.names)
    cols[-1] = 'Geology'
    ref_conc31.index.names = cols
    ref_conc31['upper_stdev'] = ref_conc31['up'] + ref_conc31['mean']
    ref_conc31['lower_stdev'] = ref_conc31['mean'] - ref_conc31['down']

    # Convert to log scale
    ref_conc31['log_mean'] = np.log(ref_conc31['mean'])
    ref_conc31['log_stdev'] = np.log(ref_conc31['upper_stdev']) - np.log(ref_conc31['mean'])

    ref_conc32 = ref_conc31[['log_mean', 'log_stdev']].stack()
    cols = list(ref_conc32.index.names)
    cols[-1] = 'stat'
    ref_conc32.index.names = cols

    ref_conc33 = ref_conc32.unstack(0).reset_index()

    ref_conc00 = pd.merge(seg_class0, ref_conc33, on=['Climate', 'Topography', 'Geology'], how='left')
    # ref_conc000 = ref_conc00.loc[~ref_conc00['DRP'].isnull()]

    ref_conc20 = pd.read_csv(params.rivers_ref_conc2_csv_path)
    ref_conc20['Climate'] = ref_conc20.REC.str[:2]
    ref_conc20['Topography'] = ref_conc20.REC.str[2:]

    ref_conc21 = ref_conc20.set_index(['param', 'Climate', 'Topography']).drop('REC', axis=1)
    ref_conc21['upper_stdev'] = ref_conc21['up'] + ref_conc21['mean']
    ref_conc21['lower_stdev'] = ref_conc21['mean'] - ref_conc21['down']

    # Convert to log scale
    ref_conc21['log_mean'] = np.log(ref_conc21['mean'])
    ref_conc21['log_stdev'] = np.log(ref_conc21['upper_stdev']) - np.log(ref_conc21['mean'])

    ref_conc22 = ref_conc21[['log_mean', 'log_stdev']].stack()
    cols = list(ref_conc22.index.names)
    cols[-1] = 'stat'
    ref_conc22.index.names = cols

    ref_conc23 = ref_conc22.unstack(0).reset_index()

    ref_conc000 = ref_conc00.loc[ref_conc00['DRP'].isnull(), ['nzsegment', 'Climate', 'Topography', 'Geology' ,'stat']]
    ref_conc001_list = []
    for i, grp in ref_conc000.groupby(['nzsegment', 'Climate', 'Topography', 'Geology']):
        for stat in ['log_mean', 'log_stdev']:
            l1 = list(i)
            l1.append(stat)
            ref_conc001_list.append(l1)

    ref_conc001 = pd.DataFrame(ref_conc001_list, columns=['nzsegment', 'Climate', 'Topography', 'Geology', 'stat'])

    ref_conc01 = pd.merge(ref_conc001, ref_conc23, on=['Climate', 'Topography', 'stat'], how='left')

    ref_conc24 = ref_conc23.groupby(['Climate', 'stat']).mean()
    other1 = ref_conc24.groupby('stat').mean().reset_index()
    other1['Climate'] = 'Other'
    ref_conc25 = pd.concat([ref_conc24, other1.set_index(['Climate', 'stat'])])

    ref_conc02 = pd.merge(ref_conc01.loc[ref_conc01['DRP'].isnull(), ['nzsegment', 'Climate', 'Topography', 'Geology', 'stat']], ref_conc25, on=['Climate', 'stat'], how='left')

    ref_conc000 = ref_conc00.loc[~ref_conc00['DRP'].isnull()]
    ref_conc010 = ref_conc01.loc[~ref_conc01['DRP'].isnull()]

    ref_conc0 = pd.concat([ref_conc000, ref_conc010, ref_conc02]).drop(['Climate', 'Topography', 'Geology'], axis=1).set_index(['nzsegment', 'stat']).round(4).unstack(1)

    ref_conc0.to_csv(params.rivers_ref_conc_csv_path)

    ref_conc0.columns.names = ['indicator', 'stat']

    ## Renaming parameters
    # ref_conc0 = ref_conc0.rename(columns={'NO3N': 'NNN', 'SS': 'sediment', 'ECOLI': 'e.coli'})

    # cols2 = ref_conc0.columns
    # for param, col in utils.indicator_dict.items():
    #     ref_conc0[param] = ref_conc0[col].copy()

    # ref_conc0['Ammoniacal nitrogen'] = ref_conc0['NH4N'].copy()

    # ref_conc1 = ref_conc0.drop(cols2, axis=1)

    # ref_conc1 = ref_conc0

    # ref_load0 = pd.merge(flows, ref_conc1, on='nzsegment')
    # cols1 = [col for col in ref_load0.columns if col not in ['nzsegment', 'flow']]
    # for col in cols1:
    #     ref_load0[col] = ref_load0[col] * ref_load0['flow']

    # ref_load1 = ref_load0.drop(['flow'], axis=1).set_index('nzsegment')
    # ref_load1.columns = [col+'_reference_load' for col in ref_load1.columns]
    # ref_load1.to_csv(params.rivers_ref_load_csv_path)

    ## Save as dict
    # ref_conc_dict = {}
    # for i, row in ref_conc0.unstack(1).iterrows():
    #     ref_conc_dict[i] = row.unstack(0).to_dict()

    # ref_conc_dict = ref_conc0.to_dict('index')

    ## Yields per segment
    q_mean0 = pd.read_csv(params.rivers_loads_csv, usecols=['nzsegment', 'Qmean']).set_index('nzsegment').rename(columns={'Qmean': 'q_mean'})['q_mean']

    q_segs = set(q_mean0.index)

    q_mean1 = q_mean0.to_dict()

    with booklet.open(params.rec_reference_conc_path, 'n', value_serializer='pickle', key_serializer='uint4', n_buckets=10007) as f:
        with booklet.open(params.lake_outflows_blt) as out:
            for LFENZID, segs in out.items():
                # if len(segs) > 1:
                #     d
                p_tot = {}
                q_tot = 0
                for seg in segs:
                    if seg in q_segs:
                        q = q_mean1[seg]
                        q_tot += q

                        ref0 = ref_conc0.loc[seg]

                        for ind, val in ref0.groupby('indicator'):
                            val0 = val.loc[ind]
                            mean1 = val0['log_mean']
                            stdev = val0['log_stdev']
                            ref_5 = mean1 - (stdev * 1.65)
                            ref_95 = mean1 + (stdev * 1.65)
                            ref_5_conc = np.exp(ref_5) * q
                            ref_95_conc = np.exp(ref_95) * q
    
                            if ind in p_tot:
                                p_tot[ind][5] += ref_5_conc
                                p_tot[ind][95] += ref_95_conc
                                p_tot[ind][50] += np.exp(mean1) * q
                            else:
                                p_tot[ind] = {5: ref_5_conc, 95: ref_95_conc, 50: np.exp(mean1) * q}

                if p_tot:
                    res = {}
                    for ind, data in p_tot.items():
                        res[ind] = {5: round(data[5]/q_tot, 4), 95: round(data[95]/q_tot, 4), 50: round(data[50]/q_tot, 4)}
                    f[LFENZID] = res

    ### Lake reference values
    ref_conc0 = pd.read_csv(params.lake_ref_conc_csv).dropna()

    # Missing lakes
    lakes_poly1 = gpd.read_file(params.lakes_poly_gpkg_path)
    ref_conc1 = ref_conc0[ref_conc0.LFENZID.isin(lakes_poly1.LFENZID)].copy()
    lakes_poly_missing = lakes_poly1[~lakes_poly1.LFENZID.isin(ref_conc1.LFENZID)].copy()
    lakes_poly_missing['area_ha'] = (lakes_poly_missing.area/10000).round(1)
    lakes_poly_missing.drop('geometry', axis=1).to_csv(params.lake_missing_ref_csv, index=False)

    ## process data
    ref_conc1 = ref_conc0.drop('Region', axis=1).set_index('LFENZID')
    ref_conc1.columns = [c.split('reference_')[1] for c in ref_conc1.columns]
    ref_conc1 = ref_conc1.round(3).rename(columns={'SECCHI': 'Secchi'})

    with booklet.open(params.lake_reference_conc_path, 'n', value_serializer='orjson', key_serializer='uint4', n_buckets=10007) as f:
        for lake_id, val in ref_conc1.to_dict('index').items():
            f[lake_id] = val

























































