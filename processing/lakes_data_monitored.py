#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 14:53:57 2023

@author: mike
"""
import os
import pandas as pd
import numpy as np
from statsmodels.tsa import seasonal
import scipy
import geopandas as gpd
import booklet

import utils, params


pd.options.display.max_columns = 10


#############################################################
### Parameters

# params = ['TN', 'Secchi', 'TP', 'Chla']
# base_dir_name = '{param}_Decomposed_Timeseries'
# base_file_name = 'Decompose'
# freq_code = 'Q'

#############################################################
### Functions


def reg_transform(x, y, slope, intercept):
    """

    """
    y_new = y - (slope*x + intercept)
    return y_new


def est_median_freq(data):
    """

    """
    d2 = data['date'].shift(1) - data['date']
    m_days = -d2.dt.days.median()

    return m_days


def regress(results, data_col, date_col='date'):
    """

    """
    grp2 = results.groupby(['parameter', 'lawa_id'])

    r_list = []
    for i, data in grp2:
        x = data[date_col].values.astype('datetime64[D]').astype(int)
        y = data[data_col].values
        slope, intercept, r, p, se = scipy.stats.linregress(x, y)
        if p < 0.05:
            new_y = []
            for x1, y1 in zip(x, y):
                new_y.append(reg_transform(x1, y1, slope, intercept))
            new_y = np.array(new_y)
        else:
            new_y = y - y.mean()

        df1 = pd.DataFrame(zip(data[date_col].values, new_y), columns=['date', data_col])
        df1['parameter'] = i[0]
        df1['lawa_id'] = i[1]
        r_list.append(df1)

    df2 = pd.concat(r_list).set_index(['parameter', 'lawa_id', 'date'])

    return df2


def deseason(data):
    """

    """
    freq_code = '2M'
    grp1 = data.groupby(['lawa_id', 'parameter'])[['date', 'observed']]

    raw_output_list = []
    for i, data in grp1:
        d1 = data.set_index('date')['observed']
        # d1 = d1[d1 != 0]
        # d2 = d1.reset_index()['date'].shift(1) - d1.index
        # m_days = -d2.dt.days.median()
        # if m_days < 60:
        #     freq_code = 'M'
        #     # seasonal = 13
        # elif m_days < 90:
        #     freq_code = '2M'
        #     # seasonal = 7
        # else:
        #     freq_code = 'Q'
        reg1 = pd.date_range(d1.index[0], d1.index[-1], freq=freq_code)
        reg2 = reg1[~reg1.isin(d1.index)]
        s1 = pd.Series(np.nan, index=reg2)
        s2 = pd.concat([d1, s1]).sort_index()
        # s2 = pd.concat([np.log(d1), s1]).sort_index()
        s3 = s2.interpolate('time')
        s4 = (s3 + s3.shift(-1))/2
        s5 = s4.resample(freq_code).mean().dropna()
        # s5 = s3[reg1]
        s5.name = 'observed'

        r1 = seasonal.STL(s5, robust=False, seasonal=13).fit()
        r2 = pd.concat([r1.observed, r1.trend, r1.seasonal, r1.resid], axis=1)
        r2.index.name = 'date'

        # Put the observations back in and use a linear interp to get the others
        # r2b = pd.concat([np.log(d1).to_frame(), r2]).sort_index()
        r2b = pd.concat([d1.to_frame(), r2]).sort_index()
        r2b = r2b[~r2b.index.duplicated(keep='last')].copy()
        r2b[['trend', 'season']] = r2b[['trend', 'season']].interpolate('time', limit_area='inside')
        r2b['resid'] = r2b['observed'] - r2b['trend'] - r2b['season']
        r2d = r2b.loc[d1.index, :].dropna()

        r2d['lawa_id'] = i[0]
        r2d['parameter'] = i[1]
        r3 = r2d.reset_index().set_index(['parameter', 'lawa_id', 'date'])
        raw_output_list.append(r3)

    all_results = pd.concat(raw_output_list)

    return all_results


def test_deseason_resampling(data):
    """

    """
    grp1 = data.groupby(['lawa_id', 'parameter'])[['date', 'observed']]

    output_list = []
    for i, data in grp1:
        d1 = data.set_index('date')['observed']
        d2 = d1.reset_index()['date'].shift(1) - d1.index
        m_days = -d2.dt.days.median()
        if m_days < 90:
            train = d1[::2]
            test = d1[1::2]
            d2 = train.reset_index()['date'].shift(1) - train.index
            m_days = -d2.dt.days.median()
            if m_days < 60:
                freq_code = '1M'
            elif m_days < 90:
                freq_code = '2M'
            else:
                freq_code = 'Q'

            ## Use full dataset at lower freq
            reg1 = pd.date_range(d1.index[0], d1.index[-1], freq=freq_code)
            reg2 = reg1[~reg1.isin(d1.index)]
            s1 = pd.Series(np.nan, index=reg2)
            s2 = pd.concat([d1, s1]).sort_index()
            s3 = s2.interpolate('time')
            s4 = (s3 + s3.shift(-1))/2
            s5 = s4.resample(freq_code).mean().dropna()
            s5.name = 'observed'

            r1 = seasonal.STL(s5, robust=False, seasonal=13).fit()
            r2 = pd.concat([r1.observed, r1.trend, r1.seasonal, r1.resid], axis=1)
            r2.index.name = 'date'

            # Put the observations back in and use a linear interp to get the others
            r2b = pd.concat([test.to_frame(), r2]).sort_index()
            r2b = r2b[~r2b.index.duplicated(keep='last')].copy()
            r2b[['trend', 'season']] = r2b[['trend', 'season']].interpolate('time', limit_area='inside')
            r2b['resid'] = r2b['observed'] - r2b['trend'] - r2b['season']
            r2d = r2b.loc[test.index, :].dropna()
            full1 = r2d['trend'] + r2d['resid']
            # base_stdev = r2e.std()

            ## Test by removing half the data
            reg1 = pd.date_range(train.index[0], train.index[-1], freq=freq_code)
            reg2 = reg1[~reg1.isin(train.index)]
            s1 = pd.Series(np.nan, index=reg2)
            s2 = pd.concat([train, s1]).sort_index()
            s3 = s2.interpolate('time')
            s4 = (s3 + s3.shift(-1))/2
            s5 = s4.resample(freq_code).mean().dropna()
            s5.name = 'observed'

            r1 = seasonal.STL(s5, robust=False, seasonal=13).fit()
            r2 = pd.concat([r1.observed, r1.trend, r1.seasonal, r1.resid], axis=1)
            r2.index.name = 'date'

            # Put the observations back in and use a linear interp to get the others
            r2b = pd.concat([test.to_frame(), r2]).sort_index()
            r2b = r2b[~r2b.index.duplicated(keep='last')].copy()
            r2b[['trend', 'season']] = r2b[['trend', 'season']].interpolate('time', limit_area='inside')
            r2b['resid'] = r2b['observed'] - r2b['trend'] - r2b['season']
            r2d = r2b.loc[test.index, :].dropna()
            test1 = r2d['trend'] + r2d['resid']
            # test_stdev = r2e.std()

            ## Combine results
            combo3 = pd.concat([full1, test1], axis=1).dropna()
            base_stdev, test_stdev = combo3.std()
            r3 = [i[1], i[0], freq_code, base_stdev, test_stdev]
            output_list.append(r3)

    all_results = pd.DataFrame(output_list, columns=['parameter', 'lawa_id', 'freq', 'base_stdev', 'test_stdev'])
    all_results['error'] = (all_results['test_stdev'] / all_results['base_stdev']) - 1
    summ_results1 = all_results.groupby(['parameter', 'freq']).mean(numeric_only=True)
    summ_results2 = all_results.groupby(['parameter', 'freq'])['error'].count()
    summ_results2.name = 'site_count'
    summ_results = pd.concat([summ_results1, summ_results2], axis=1)

    return all_results.set_index(['parameter', 'lawa_id', 'freq']), summ_results


def test_deseason_resampling_2M(data):
    """

    """
    freq_code = '2M'
    grp1 = data.groupby(['lawa_id', 'parameter'])[['date', 'observed']]

    output_list = []
    for i, data1 in grp1:
        d1 = data1.set_index('date')['observed']

        base_stdev_sum = pd.Series([0]*5, index=['observed', 'trend', 'season', 'resid', 'trend+resid'])
        test_stdev_sum = base_stdev_sum.copy()
        for t in range(2):
            if t == 0:
                train_start = 0
                test_start = 1
            else:
                train_start = 1
                test_start = 0
            train = d1[train_start::2]
            test = d1[test_start::2]
            # d2 = d1.reset_index()['date'].shift(1) - d1.index
            # m_days = -d2.dt.days.median()
            # if m_days < 90:
            #     train = d1[::2]
            #     test = d1[1::2]
            #     d2 = train.reset_index()['date'].shift(1) - train.index
            #     m_days = -d2.dt.days.median()
            #     if m_days < 60:
            #         freq_code = '1M'
            #     elif m_days < 90:
            #         freq_code = '2M'
            #     else:
            #         freq_code = 'Q'

            ## Use full dataset at lower freq
            reg1 = pd.date_range(d1.index[0], d1.index[-1], freq=freq_code)
            reg2 = reg1[~reg1.isin(d1.index)]
            s1 = pd.Series(np.nan, index=reg2)
            s2 = pd.concat([d1, s1]).sort_index()
            s3 = s2.interpolate('time')
            s4 = (s3 + s3.shift(-1))/2
            s5 = s4.resample(freq_code).mean().dropna()
            s5.name = 'observed'

            r1 = seasonal.STL(s5, robust=False, seasonal=13).fit()
            r2 = pd.concat([r1.observed, r1.trend, r1.seasonal, r1.resid], axis=1)
            r2.index.name = 'date'

            # Put the observations back in and use a linear interp to get the others
            r2b = pd.concat([test.to_frame(), r2]).sort_index()
            r2b = r2b[~r2b.index.duplicated(keep='last')].copy()
            r2b[['trend', 'season']] = r2b[['trend', 'season']].interpolate('time', limit_area='inside')
            r2b['resid'] = r2b['observed'] - r2b['trend'] - r2b['season']
            r2d = r2b.loc[test.index, :].dropna()
            r2d['trend+resid'] = r2d['trend'] + r2d['resid']
            # full1 = r2d['trend'] + r2d['resid']
            base_stdev = r2d.std()

            ## Test by removing half the data
            reg1 = pd.date_range(train.index[0], train.index[-1], freq=freq_code)
            reg2 = reg1[~reg1.isin(train.index)]
            s1 = pd.Series(np.nan, index=reg2)
            s2 = pd.concat([train, s1]).sort_index()
            s3 = s2.interpolate('time')
            s4 = (s3 + s3.shift(-1))/2
            s5 = s4.resample(freq_code).mean().dropna()
            s5.name = 'observed'

            r1 = seasonal.STL(s5, robust=False, seasonal=13).fit()
            r2 = pd.concat([r1.observed, r1.trend, r1.seasonal, r1.resid], axis=1)
            r2.index.name = 'date'

            # Put the observations back in and use a linear interp to get the others
            r2b = pd.concat([test.to_frame(), r2]).sort_index()
            r2b = r2b[~r2b.index.duplicated(keep='last')].copy()
            r2b[['trend', 'season']] = r2b[['trend', 'season']].interpolate('time', limit_area='inside')
            r2b['resid'] = r2b['observed'] - r2b['trend'] - r2b['season']
            r2d = r2b.loc[test.index, :].dropna()
            r2d['trend+resid'] = r2d['trend'] + r2d['resid']
            # test1 = r2d['trend'] + r2d['resid']
            test_stdev = r2d.std()

            ## Combine results
            # combo3 = pd.concat([full1, test1], axis=1).dropna()
            # base_stdev, test_stdev = combo3.std()
            base_stdev_sum += base_stdev
            test_stdev_sum += test_stdev

        base_stdev_sum1 = base_stdev_sum/2
        # base_stdev_sum1['type'] = 'base'
        # base_stdev_sum1['parameter'] = i[1]
        # base_stdev_sum1['lawa_id'] = i[0]
        test_stdev_sum1 = test_stdev_sum/2
        # test_stdev_sum1['type'] = 'test'
        # test_stdev_sum1['parameter'] = i[1]
        # test_stdev_sum1['lawa_id'] = i[0]

        stdev_sum1 = (test_stdev_sum1/base_stdev_sum1) - 1
        stdev_sum1['parameter'] = i[1]
        stdev_sum1['lawa_id'] = i[0]

        # r3 = [i[1], i[0], base_stdev_sum/2, test_stdev_sum/2]
        output_list.append(stdev_sum1.to_frame().transpose())

    # all_results = pd.DataFrame(output_list, columns=['parameter', 'lawa_id', 'base_stdev', 'test_stdev'])
    all_results = pd.concat(output_list).reset_index(drop=True)
    for col in ['observed', 'trend', 'season', 'resid', 'trend+resid']:
        all_results[col] = all_results[col].astype(float)
    # all_results['error'] = (all_results['test_stdev'] / all_results['base_stdev']) - 1
    summ_results1 = all_results.groupby(['parameter']).mean(numeric_only=True)
    summ_results2 = all_results.groupby(['parameter'])['lawa_id'].count()
    summ_results2.name = 'site_count'
    summ_results = pd.concat([summ_results1, summ_results2], axis=1)

    return all_results.set_index(['parameter', 'lawa_id']), summ_results


##############################################################
### Process data

## monitoring data from individual files
# data_list = []
# for param in params:
#     dir_name = base_dir_name.format(param=param)
#     param_path = utils.lakes_source_path.joinpath(dir_name)
#     for path in param_path.iterdir():
#         # site_name = path.name.split(base_file_name)[1].split('.')[0]
#         data = pd.read_csv(path).iloc[:, 1:].rename(columns={'ID': 'site_id', 'Date': 'date', 'Observed': 'observed'})
#         data['date'] = pd.to_datetime(data['date'], dayfirst=True)
#         data['parameter'] = param
#         data_list.append(data)

# data0 = pd.concat(data_list).set_index(['site_id', 'parameter', 'date']).sort_index()
# data0.to_csv(utils.lakes_source_data_path)


def lakes_monitored_conc():
    ## monitoring data from large spreadsheet
    ts_data0 = pd.read_feather(params.lake_moni_data_raw)
    site_data0 = gpd.read_file(params.lake_site_data_raw)
    mtypes0 = pd.read_csv(params.lake_mtypes_data_raw)

    site_data1 = site_data0[site_data0.LFENZID > 0].copy()
    site_data1['LFENZID'] = site_data1['LFENZID'].astype('int32')
    ts_data1 = ts_data0[ts_data0.LawaSiteID.isin(site_data1.LawaSiteID)].copy()
    ts_data1.loc[ts_data1.Symbol == 'Right', 'Symbol'] = 'greater_than'
    ts_data1.loc[ts_data1.Symbol == 'Left', 'Symbol'] = 'less_than'
    ts_data1.loc[~ts_data1.Symbol.isin(['greater_than', 'less_than']), 'Symbol'] = 'not_censored'

    ## Filter out sites with less than 24 measurements
    mtypes1 = mtypes0[mtypes0.value_count > 24].copy()
    ts_data2 = pd.merge(mtypes1.drop('value_count', axis=1), ts_data1, on=['LawaSiteID', 'Indicator'])

    ## Only include specific indicators
    ts_data2 = ts_data2[ts_data2['Indicator'].isin(params.indicators)].copy()

    ## Convert concentrations to mg/m3
    ts_data2.loc[ts_data2.Units == 'mg/L', 'Value'] = ts_data2.loc[ts_data2.Units == 'mg/L', 'Value'] * 1000

    ## Stats on site above/below dtl
    dtl_list = []
    for i, df in ts_data2.groupby(['LawaSiteID', 'Indicator']):
        dtl_bool = df.Symbol.isin(['greater_than', 'less_than'])
        # if dtl_bool.any():
        ratio = round(dtl_bool.sum()/len(dtl_bool), 3)
        dtl_list.append([i[0], i[1], ratio])

    dtl_df0 = pd.DataFrame(dtl_list, columns=['LawaSiteID', 'Indicator', 'dtl_ratio'])

    dtl_df1 = dtl_df0[dtl_df0.dtl_ratio <= 0.4]

    ts_data3 = pd.merge(dtl_df1[['LawaSiteID', 'Indicator']], ts_data2, on=['LawaSiteID', 'Indicator']).drop(['Units', 'QCRaw', 'QCNumber', 'QCNEMSEquivalent', 'RawValue'], axis=1).rename(columns={'LawaSiteID': 'lawa_id', 'Indicator': 'indicator', 'SampleDateTime': 'date', 'Value': 'value', 'Symbol': 'censor_code'})

    ## Convert dtls
    moni5 = utils.dtl_correction(ts_data3, 'half').drop('censor_code', axis=1)

    ## Must have recent data
    moni6 = moni5.set_index(['lawa_id', 'indicator'])
    grp = moni6.groupby(['lawa_id', 'indicator'])

    max_dates = grp.date.max()
    recent_bool = (max_dates > '2021-07-01')
    # recent_bool.name = 'remove'

    max_dates = max_dates[recent_bool] + pd.DateOffset(seconds=1)
    max_dates.name = 'max_date'

    ## Must have at least 20 samples in the past 5 years
    min_dates = pd.to_datetime((max_dates - pd.DateOffset(years=5) - pd.offsets.MonthBegin()).dt.date)
    min_dates.name = 'min_date'
    min_max_dates = pd.concat([min_dates, max_dates], axis=1)

    data_list = []
    for i, data in min_max_dates.iterrows():
        data1 = moni6.loc[i]
        data2 = data1[(data1.date >= data.min_date) & (data1.date <= data.max_date)].reset_index()
        if len(data2) >= 20:
            data_list.append(data2)

    wq_data1 = pd.concat(data_list)

    ## check the sampling frequency and only use data with < 93 days freq
    freq_list = []
    for i, row in wq_data1.groupby(['lawa_id', 'indicator']):
        freq = est_median_freq(row)
        if freq < 93:
            freq_list.append([*i, freq])

    freq_df = pd.DataFrame(freq_list, columns=['lawa_id', 'indicator', 'freq'])

    wq_data2 = pd.merge(freq_df[['lawa_id', 'indicator']], wq_data1, on=['lawa_id', 'indicator'])

    ## Filter the sites
    site_data2 = site_data1[site_data1['LawaSiteID'].isin(wq_data2.lawa_id.unique())].rename(columns={'LawaSiteID': 'lawa_id'})

    ## Save results
    site_data2.to_file(params.lake_site_data)
    utils.df_to_feather(moni5, params.lake_moni_data)

    moni6 = pd.merge(site_data2[['lawa_id', 'LFENZID']], moni5, on='lawa_id')
    moni7 = moni6.groupby(['LFENZID', 'indicator', 'date'])['value'].mean()
    lake_ids = moni7.index.get_level_values(0).unique()

    # Save data for app
    with booklet.open(params.lake_moni_conc_data_blt, 'n', value_serializer='pickle', key_serializer='uint4', n_buckets=10007) as f:
        for LFENZID in lake_ids:
            f[LFENZID] = moni7.loc[LFENZID]

    ## Calculate the medians per site and lake
    site_median1 = wq_data2.groupby(['lawa_id', 'indicator'])['value'].median()

    lake_median0 = pd.merge(site_data2[['lawa_id', 'LFENZID']], site_median1.reset_index(), on='lawa_id')
    # lake_median1 = lake_median0.groupby(['LFENZID', 'indicator'])['value'].mean()
    lake_median0.to_csv(params.lake_moni_conc_csv, index=False)

    lake_median1 = lake_median0.groupby(['LFENZID', 'indicator'])['value'].mean()
    lake_ids = lake_median1.index.get_level_values(0).unique()

    # Save data for app
    with booklet.open(params.lake_moni_conc_blt, 'n', value_serializer='orjson', key_serializer='uint4', n_buckets=10007) as f:
        for LFENZID in lake_ids:
            f[LFENZID] = lake_median1.loc[LFENZID].to_dict()











































