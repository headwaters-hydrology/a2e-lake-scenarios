#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 11:01:00 2023

@author: mike
"""
import pathlib
import zstandard as zstd
import xarray as xr
import pickle
import codecs
import numpy as np
import geopandas as gpd
import base64
import io
import booklet
import pandas as pd
import plotly.graph_objects as go
from dash import dcc, html, dash_table
from copy import deepcopy

# import parameters as param
import utils.parameters as param

##############################################
### Parameters




###############################################
### Helper Functions


def sep_reference_values(ref_str):
    """

    """
    mean1, lower_upper = ref_str.split(' [ ')
    lower1, upper1 = lower_upper[:-2].split(' - ')

    return float(mean1), float(lower1), float(upper1)


def get_value(file_path, key):
    """

    """
    with booklet.open(file_path) as f:
        value = f[key]

    return value


def encode_obj(obj):
    """

    """
    cctx = zstd.ZstdCompressor(level=1)
    c_obj = codecs.encode(cctx.compress(pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL)), encoding="base64").decode()

    return c_obj


def decode_obj(str_obj):
    """

    """
    dctx = zstd.ZstdDecompressor()
    obj1 = dctx.decompress(codecs.decode(str_obj.encode(), encoding="base64"))
    d1 = pickle.loads(obj1)

    return d1


def calc_scenario_results(site_id, indicator, nzsegment, lc_tbl, calc_ready, trig):
    """

    """
    conc_factor = 1

    if indicator != 'ECOLI':
        data = get_value(param.sites_lc_red_yields_path, nzsegment).reindex(param.lu_order)
    
        if indicator in ('DRP', 'TP'):
            name = 'phosphorus_'
        elif indicator in ('NO3N', 'TN'):
            name = 'nitrogen_'
    
        # cols = ['land_cover', 'area_ha']
        # up_cols = ['land_cover', 'area_ha']
        cols = ['area_ha']
        up_cols = ['area_ha']
        for col in data.columns:
            if name in col:
                cols.append(col)
                up_cols.append(col[len(name):])
    
        data0 = data[cols].copy()
        data0.columns = up_cols
        # data0 = data0.set_index('land_cover')
    
        if (trig == 'tbl_calc_btn') and (calc_ready == 1):
            data = pd.DataFrame(lc_tbl).set_index('land_cover').astype('int16')
            data2 = pd.concat([data0, data], axis=1)
            # print(data2)
            new_contrib = (data2['yield']*data2['new_land_area'])*(1 - data2['mitigation']*0.01)
            data2['new_load'] = ((new_contrib/new_contrib.sum()) * 100).astype('int8')
            old_contrib = data2['yield']*data2['land_area']
            conc_factor = new_contrib.sum()/old_contrib.sum()
    
            old_tbl_data = deepcopy(lc_tbl)
            lc_tbl = []
            for ld in old_tbl_data:
                lc = ld['land_cover']
                ld['new_load'] = data2.loc[lc, 'new_load']
                lc_tbl.append(ld)

        else:
            tot_area = data0['area_ha'].sum()
            area_perc = ((data0['area_ha']/tot_area) * 100).round().astype('int8')
            sum1 = area_perc.sum()
            if sum1 != 100:
                diff = 100 - sum1
                area_perc.loc[area_perc.idxmax()] = area_perc.loc[area_perc.idxmax()] + diff
    
            mass_ish = (data0['yield']*area_perc)
            tot_mass = mass_ish.sum()
            mass_perc = ((mass_ish/tot_mass) * 100).round().astype('int8')
            sum1 = mass_perc.sum()
            if sum1 != 100:
                diff = 100 - sum1
                mass_perc.loc[mass_perc.idxmax()] = mass_perc.loc[mass_perc.idxmax()] + diff
            reductions = data0['reduction'].round()
    
            lc_tbl = []
            for i, lc in enumerate(data0.index):
                dict1 = {'land_cover': lc, 'land_area': area_perc[i], 'load': mass_perc[i], 'mitigation': reductions[i], 'new_land_area': area_perc[i]}
    
                lc_tbl.append(dict1)

    else:
        data = get_value(param.ecoli_catch_yields_path, nzsegment)

        if (trig == 'tbl_calc_btn') and (calc_ready == 1):
            # print(lc_tbl)
            data0 = pd.DataFrame(lc_tbl).set_index('land_cover').astype('int8') * 0.01
            data0.index.name = 'land_use'
            mitigation = data0[['mitigation']]
            # print(data0)
            # data0.loc['Dairy', 'new_land_area'] = data0.loc['Dairy', 'new_land_area'] - 0.1
            # data0.loc['Sheep and Beef', 'new_land_area'] = data0.loc['Sheep and Beef', 'new_land_area'] - 0.2
            # data0.loc['Exotic Forest', 'new_land_area'] = data0.loc['Exotic Forest', 'new_land_area'] + 0.25
            # data0.loc['Native Vegetation', 'new_land_area'] = data0.loc['Native Vegetation', 'new_land_area'] + 0.05
            lc_ratio_change = (data0['new_land_area'] - data0['land_area']).round(2)
            lc_ratio_change.name = 'area_change_ratio'
            areas_lost1 = lc_ratio_change[lc_ratio_change < 0] * -1
            if not areas_lost1.empty:
                areas_lost2 = data[data.land_use.isin(areas_lost1.index)]
                areas_lost3 = pd.merge(areas_lost2, areas_lost1, on='land_use')
                areas_lost3['area_change'] = areas_lost3['area_change_ratio'] * areas_lost3['area_ha']
                areas_lost0 = areas_lost3.groupby(['elev', 'drain'])['area_change'].sum().reset_index()
                areas_lost3['area_ha'] = areas_lost3['area_ha'] - areas_lost3['area_change']

                tot_change = areas_lost1.sum()

                areas_gained1 = lc_ratio_change[lc_ratio_change > 0]
                areas_gained2 = data[data.land_use.isin(areas_gained1.index)]
                areas_gained3 = pd.merge(areas_gained2, areas_gained1, on='land_use')
                areas_gained4 = pd.merge(areas_gained3, areas_lost0, on=['elev', 'drain'])
                areas_gained4['area_change'] = ((areas_gained4['area_change_ratio']/tot_change) * areas_gained4['area_change'])
                areas_gained4['area_ha'] = areas_gained4['area_ha'] + areas_gained4['area_change']

                areas_unchanged1 = lc_ratio_change[lc_ratio_change == 0]
                areas_unchanged2 = data[data.land_use.isin(areas_unchanged1.index)]

                combo1 = pd.concat([areas_unchanged2, areas_lost3.drop(['area_change_ratio', 'area_change'], axis=1), areas_gained4.drop(['area_change_ratio', 'area_change'], axis=1)])

                combo1 = pd.merge(combo1, mitigation, on='land_use')
                combo1['area_ratio'] = combo1['area_ha']/combo1['area_ha'].sum()
                combo1['ecoli_load'] = combo1['area_ratio'] * combo1['ecoli_factor'] * (1 - combo1['mitigation'])

                combo2 = calc_ecoli_load(combo1)

                old_tbl_data = deepcopy(lc_tbl)
                lc_tbl = []
                for ld in old_tbl_data:
                    lc = ld['land_cover']
                    ld['new_load'] = combo2.loc[lc, 'ecoli_load_perc']
                    lc_tbl.append(ld)

            else:
                combo1 = pd.merge(data, mitigation, on='land_use')
                combo1['area_ratio'] = combo1['area_ha']/combo1['area_ha'].sum()
                combo1['ecoli_load'] = combo1['area_ratio'] * combo1['ecoli_factor'] * (1 - combo1['mitigation'])
                combo1['ecoli_load_perc'] = ((combo1['ecoli_load']/combo1['ecoli_load'].sum()) * 100).astype('int8')

                combo2 = combo1.groupby('land_use')[['ecoli_load_perc']].sum()
                combo2.index.name = 'land_cover'

                old_tbl_data = deepcopy(lc_tbl)
                lc_tbl = []
                for ld in old_tbl_data:
                    lc = ld['land_cover']
                    ld['new_load'] = combo2.loc[lc, 'ecoli_load_perc']
                    lc_tbl.append(ld)

            ## Calc conc factor for indicator
            data['area_ratio'] = data['area_ha']/data['area_ha'].sum()
            data['ecoli_load'] = data['area_ratio'] * data['ecoli_factor']
            old_ecoli_load = data['ecoli_load'].sum()
            new_ecoli_load = combo1['ecoli_load'].sum()
            conc_factor = new_ecoli_load/old_ecoli_load

        else:
            data['area_ratio'] = data['area_ha']/data['area_ha'].sum()
            data['ecoli_load'] = data['area_ratio'] * data['ecoli_factor']
            combo2 = calc_ecoli_load(data)

            lc_tbl = []
            for lc, row in combo2.iterrows():
                dict1 = {'land_cover': lc, 'land_area': row['area_ratio'], 'load': row['ecoli_load_perc'], 'mitigation': int(param.ecoli_reductions[lc] * 100), 'new_land_area': row['area_ratio']}
    
                lc_tbl.append(dict1)

    return lc_tbl, conc_factor



#########################################
### Table creation


def make_results_table(data=[], tab='ind', indicator='ECOLI'):
    """

    """
    if tab == 'ind':
        if indicator == 'ECOLI':
            cols_dict = param.ind_tbl_cols_ecoli_dict
        else:
            cols_dict = param.ind_tbl_cols_main_dict

        tbl = dash_table.DataTable(
            data=data,
            style_table={'overflowX': 'auto'},
            style_cell={
                'whiteSpace': 'normal',
                'height': 'auto',
                'font-family':'sans-serif',
                'font-size': 12
                },
            columns=[{'name': name, 'id': key} for key, name in cols_dict.items()],
            # id='stats_tbl',
            style_header_conditional=[{
                'if': {'column_id': 'conc'},
                'font-weight': 'bold'
            },
                {
                'if': {
                    # 'filter_query': "{name} = Target conc",
                    'column_id': 'name'},
                'font-weight': 'bold'
            },

            ],
            )
    else:
        tbl = dash_table.DataTable(
            data=data,
            style_table={'overflowX': 'auto'},
            style_cell={
                'whiteSpace': 'normal',
                'height': 'auto',
                'font-family':'sans-serif',
                'font-size': 12
                },
            columns=[{'name': name, 'id': key} for key, name in param.peri_tbl_cols_dict.items()],
            merge_duplicate_headers=True,
            # id='stats_tbl',
            style_header_conditional=[
                {
                'if': {'column_id': 'conc_shaded'},
                'font-weight': 'bold'
                },
                {
                'if': {'column_id': 'conc_unshaded'},
                'font-weight': 'bold'
                },
                {
                'if': {
                    # 'filter_query': "{name} = Target conc",
                    'column_id': 'name'},
                'font-weight': 'bold'
            },
            ],
            )

    return tbl


#############################################
### Figure creation


def make_fig_ind(stats_tbl=None, site_id=None, indicator=None):
    """
    I need to make a figure this way because of a bug in the reset axes of the figure.
    """
    stats = pd.DataFrame(stats_tbl).set_index('name')['conc']
    names = stats.index

    # print(stats)
    data = get_value(param.rivers_data_app_path, (site_id, indicator))
    ymax = data.max()
    if indicator == 'ECOLI':
        y_max_plot = np.percentile(data.values, 75)
    else:
        y_max_plot = np.percentile(data.values, 95)
    # data_len = len(data)
    start_date = data.index[0]
    end_date = data.index[-1]
    diff_days = (end_date - start_date).days
    extra_x_days = int(diff_days*0.04)
    start_x = start_date - pd.DateOffset(days=extra_x_days)
    end_x = end_date + pd.DateOffset(days=extra_x_days)

    # print(stats)

    # data_list = list(stats.values())
    # data_list.extend([stats['median']]*10)

    ## Fig
    fig = go.Figure()

    ## Add traces
    # Samples
    fig.add_trace(
        go.Scatter(
            x=data.index,
            y=data.values,
            mode='markers',
            marker={
                'color': 'rgb(179, 179, 179)',
                },
            hoverinfo='name+y',
            name='Observation',
            # legendrank=1,
            )
        )

    # Bands
    if 'Band A' in names:
        y_previous = 0
        for name in names:
            if 'Band' in name:
                if 'D' in name:
                    fig.add_trace(
                        go.Scatter(
                            x=[start_x, end_x],
                            y=[ymax - y_previous] * 2,
                            mode='lines',
                            line=dict(width=0.5, color=param.plot_colors[name]),
                            hoverinfo='name',
                            name=name,
                            stackgroup='Bands'
                            )
                        )
                else:
                    ymax1 = float(stats[name].split('- ')[1])
                    ymax_plot = ymax1 - y_previous
                    y_previous = ymax1
                    # print(ymax_plot)
                    fig.add_trace(
                        go.Scatter(
                            x=[start_x, end_x],
                            y=[ymax_plot] * 2,
                            mode='lines',
                            line=dict(width=0.5, color=param.plot_colors[name]),
                            hoverinfo='name',
                            name=name,
                            stackgroup='Bands'
                            )
                        )

    # Reference
    mean1, lower1, upper1 = sep_reference_values(stats['Reference'])

    fig.add_trace(
        go.Scatter(
            x=[start_x, end_x, end_x, start_x],
            y=[upper1, upper1, lower1, lower1],
            mode='lines',
            line=dict(
                dash='dash',
                width=0,
                color='rgb(170, 170, 170)'
                ),
            hoverinfo='name+y',
            name='Reference',
            showlegend=True,
            fill='toself',
            fillpattern=dict(
                size=10,
                solidity=0.2,
                shape='/',
                bgcolor='rgba(0, 0, 0, 0)',
                )
            )
        )

    # Current
    median1 = float(stats['Current'])
    fig.add_trace(
        go.Scatter(
            x=[start_x, end_x],
            y=[median1] * 2,
            mode='lines',
            line=dict(
                dash='dash',
                width=2,
                color=param.plot_colors['Current']
                ),
            hoverinfo='name+y',
            name='Current',
            )
        )

    # Scenario
    median1 = stats['Scenario']
    if median1 != 'Press the Run catchment scenario button':
        fig.add_trace(
            go.Scatter(
                x=[start_x, end_x],
                y=[float(median1)] * 2,
                mode='lines',
                line=dict(
                    dash='dashdot',
                    width=2,
                    color=param.plot_colors['Scenario']
                    ),
                hoverinfo='name+y',
                name='Scenario',
                )
            )

    ## Update layout
    if indicator == 'ECOLI':
        units = 'CFU/100ml'
    else:
        units = 'mg/l'
    fig.update_layout(
        yaxis_title=units,
        margin=dict(l=20, r=20, t=0, b=60),
        legend=dict(
            bgcolor="rgb(220, 220, 220)",
            bordercolor="Black",
            borderwidth=1.5,
            y=0.9,
            )
        )
    fig.update_yaxes(range=[0, y_max_plot])
    fig.update_xaxes(range=[start_x, end_x])

    return fig


def make_graph(fig=None):
    """
    I need to make a figure this way because of a bug in the reset axes of the figure.
    """
    if fig is not None:
        return html.Div(dcc.Graph(id='ts_plot', figure=fig, config={'responsive': True, 'displaylogo': False, 'modeBarButtonsToRemove': ['zoomIn', 'zoomOut', 'select2d', 'lasso2d']}, style={'height': '50vh'}), id='ts_plot_div', style={'width': '100%', 'height': 'auto', 'margin': "auto", "display": "block"})
    else:
        return html.Div(id='ts_plot_div', style={'width': '100%', 'height': 'auto', 'margin': "auto", "display": "block"})


def calc_lake_conc_change_ratio(inflow_ratio, indicator, max_depth, residence_time, ref_cond=False):
    """

    """
    if indicator == 'TP':
        if ref_cond:
            if max_depth > 7.5:
                b = 1 + 0.27*(residence_time**0.29)
                r_lake = inflow_ratio**(1/b)
            else:
                r_lake = inflow_ratio
        else:
            if max_depth > 7.5:
                b = 1 + 0.44*(residence_time**0.13)
                r_lake = inflow_ratio**(1/b)
            else:
                r_lake = inflow_ratio
    elif indicator == 'TN':
        if ref_cond:
            r_lake = inflow_ratio**0.81
        else:
            r_lake = inflow_ratio**0.54
    elif indicator == 'CHLA':
        if ref_cond:
            r_lake = inflow_ratio**1.24
        else:
            r_lake = inflow_ratio**1.25
    elif indicator == 'SECCHI':
        if ref_cond:
            r_lake = inflow_ratio**1.46
        else:
            if max_depth > 20:
                r_lake = inflow_ratio**0.9
            else:
                r_lake = inflow_ratio**0.38
    else:
        raise ValueError('No calc for indicator')

    return r_lake







































































