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

    # print(conc_factors)

    ## Nitrogen
    r_n_in = conc_factors['nitrogen'] + 0.0001
    c_n_lake_current = current_concs['TN']
    if region == 'Waikato':
        r_n_lake = est_r_lake_n_waikato(r_n_in, t)
    else:
        r_n_lake = est_r_lake_n(r_n_in)

    c_n_lake_scenario = c_n_lake_current * r_n_lake

    ## TP
    r_p_in = conc_factors['phosphorus'] + 0.0001
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
        value = f.get(key)

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


def calc_scenario_tbl(lake_id, lc_tbl, calc_ready, trig):
    """

    """
    conc_factors = None

    data = get_value(param.lakes_lc_red_yields_path, lake_id).reindex(param.lu_order)

    area_ha = data['area_ha']['area_ha']
    tot_area = area_ha.sum()

    cols = ['phosphorus', 'nitrogen']
    # if indicator in ('TP', 'CHLA', 'SECCHI'):
    #     cols.append('phosphorus')
    # elif indicator in ('TN', 'CHLA', 'SECCHI'):
    #     cols.append('nitrogen')

    # cols = ['land_cover', 'area_ha']
    # up_cols = ['land_cover', 'area_ha']
    # cols = ['area_ha']
    # up_cols = ['area_ha']
    # for col in data.columns:
    #     if name in col:
    #         cols.append(col)
    #         up_cols.append(col[len(name):])

    # data0 = data[cols].copy()
    # data0.columns = up_cols
    # data0 = data0.set_index('land_cover')

    if (trig == 'tbl_calc_btn') and (calc_ready == 1):
        data0 = pd.DataFrame(lc_tbl).set_index('land_cover').astype('int16')
        mitigation0 = data0[cols]
        # data2 = pd.concat([data0, data], axis=1)
        # print(data2)
        # new_contrib = (data['yield'][cols].multiply(tot_area * data0['new_land_area'] * 0.01, axis=0)*(1 - mitigation0*0.01))
        new_contrib = (data['yield'][cols].multiply(data0['new_land_area'], axis=0)*(1 - mitigation0*0.01))
        # data0['new_load'] = ((new_contrib.sum(axis=1)/new_contrib.sum(axis=1).sum()) * 100).astype('int8')
        old_contrib = data['yield'][cols].multiply(tot_area * data0['land_area'] * 0.01, axis=0)
        old_contrib = data['yield'][cols].multiply(data0['land_area'], axis=0)
        conc_factors = (new_contrib.sum()/old_contrib.sum()).to_dict()

        # old_tbl_data = deepcopy(lc_tbl)
        # lc_tbl = []
        # for ld in old_tbl_data:
        #     lc = ld['land_cover']
        #     ld['new_load'] = data0.loc[lc, 'new_load']
        #     lc_tbl.append(ld)

    else:
        area_perc = ((area_ha/tot_area) * 100).round().astype('int8')
        sum1 = area_perc.sum()
        if sum1 != 100:
            diff = 100 - sum1
            area_perc.loc[area_perc.idxmax()] = area_perc.loc[area_perc.idxmax()] + diff

        ## Load contributions
        mass_ish = data['yield'][cols].multiply(area_perc, axis=0)
        tot_mass = mass_ish.sum()
        mass_perc = ((mass_ish/tot_mass) * 100).round().astype('int8')
        sum1 = mass_perc.sum()
        if (sum1 != 100).any():
            diff = 100 - sum1
            mass_perc.loc[mass_perc.idxmax()] = mass_perc.loc[mass_perc.idxmax()] + diff

        ## Reductions
        reductions = data['reduction'].round().astype('int8')

        lc_tbl = []
        for i, lc in enumerate(data.index):
            red1 = reductions.iloc[i]
            dict1 = {'land_cover': lc, 'land_area': area_perc[i], 'load_n': mass_perc['nitrogen'][i], 'load_p': mass_perc['phosphorus'][i], 'nitrogen': red1['nitrogen'], 'phosphorus': red1['phosphorus'], 'new_land_area': area_perc[i]}

            lc_tbl.append(dict1)

    return lc_tbl, conc_factors



#########################################
### Table creation


def make_results_table(data=[], tab='TN'):
    """

    """
    if tab == 'Secchi':
        cols_dict = param.ind_tbl_cols_secchi_dict
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

    return tbl


#############################################
### Figure creation


def make_fig_ind(data_tbl, lake_id, indicator):
    """
    I need to make a figure this way because of a bug in the reset axes of the figure.
    """
    # stats0 = pd.DataFrame(data_tbl).set_index('name')['conc']
    # names = stats0.index

    ## Achievable concs
    achieve_stats0 = get_value(param.lake_achievable_conc_path, lake_id)
    worst = achieve_stats0['worst'][indicator]
    best = achieve_stats0['best'][indicator]

    if indicator == 'Secchi':
        y_max = best * 1.5
    else:
        y_max = worst * 1.5

    ## Fig
    fig = go.Figure()

    try:
        ## Monitoring data
        # print(stats0)
        data = get_value(param.lakes_moni_data_path, lake_id).loc[indicator]

        y_max_data = data.max()
        if indicator == 'Secchi':
            if y_max_data > y_max:
                y_max = y_max_data
        else:
            y_95th = np.percentile(data.values, 95)
            if y_95th > y_max:
                y_max = y_95th

        start_date = data.index[0]
        end_date = data.index[-1]
        diff_days = (end_date - start_date).days
        extra_x_days = int(diff_days*0.04)
        start_x = start_date - pd.DateOffset(days=extra_x_days)
        end_x = end_date + pd.DateOffset(days=extra_x_days)

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

        y_median_current = float(data_tbl['Current (measured)'])
        if data_tbl['Scenario'] == 'Press the Run catchment scenario button':
            y_median_scenario = None
        else:
            y_median_scenario = float(data_tbl['Scenario'])

        measured = True

    except:
        y_median_current = float(data_tbl['Current (modelled)'])
        if data_tbl['Scenario'] == 'Press the Run catchment scenario button':
            y_median_scenario = None
        else:
            y_median_scenario = float(data_tbl['Scenario'])

        y_max_data = y_max

        start_x = 0
        end_x = 1

        measured = False

    # print(stats)

    # data_list = list(stats.values())
    # data_list.extend([stats['median']]*10)

    # Bands
    if 'Band A' in data_tbl:
        y_previous = 0
        for name in data_tbl:
            if 'Band' in name:
                if 'D' in name:
                    fig.add_trace(
                        go.Scatter(
                            x=[start_x, end_x],
                            y=[y_max_data - y_previous] * 2,
                            mode='lines',
                            line=dict(width=0.5, color=param.plot_colors[name]),
                            hoverinfo='name',
                            name=name,
                            stackgroup='Bands'
                            )
                        )
                else:
                    ymax1 = float(data_tbl[name].split('- ')[1])
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

    # # Reference
    # ref_median = float(stats0['Reference'])

    # fig.add_trace(
    #     go.Scatter(
    #         x=[start_x, end_x],
    #         y=[ref_median] * 2,
    #         mode='lines',
    #         line=dict(
    #             dash='dash',
    #             width=8,
    #             color='rgb(170, 170, 170)'
    #             ),
    #         hoverinfo='name+y',
    #         name='Reference',
    #         showlegend=True,
    #         # fill='toself',
    #         # fillpattern=dict(
    #         #     size=10,
    #         #     solidity=0.2,
    #         #     shape='/',
    #         #     bgcolor='rgba(0, 0, 0, 0)',
    #         #     )
    #         )
    #     )

    # Worst
    # fig.add_trace(
    #     go.Scatter(
    #         x=[start_x, end_x],
    #         y=[worst] * 2,
    #         mode='lines',
    #         line=dict(
    #             dash='dash',
    #             width=4,
    #             color=param.plot_colors['worst']
    #             ),
    #         hoverinfo='name+y',
    #         name='All Dairy',
    #         showlegend=True,
    #         # fill='toself',
    #         # fillpattern=dict(
    #         #     size=10,
    #         #     solidity=0.2,
    #         #     shape='/',
    #         #     bgcolor='rgba(0, 0, 0, 0)',
    #         #     )
    #         )
    #     )

    # Best
    fig.add_trace(
        go.Scatter(
            x=[start_x, end_x],
            y=[best] * 2,
            mode='lines',
            line=dict(
                dash='dash',
                width=4,
                color=param.plot_colors['best']
                ),
            hoverinfo='name+y',
            name='All Native',
            showlegend=True,
            # fill='toself',
            # fillpattern=dict(
            #     size=10,
            #     solidity=0.2,
            #     shape='/',
            #     bgcolor='rgba(0, 0, 0, 0)',
            #     )
            )
        )

    # Current
    # median1 = float(stats0['Current'])
    fig.add_trace(
        go.Scatter(
            x=[start_x, end_x],
            y=[y_median_current] * 2,
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
    # median1 = stats0['Scenario']
    if y_median_scenario:
        fig.add_trace(
            go.Scatter(
                x=[start_x, end_x],
                y=[y_median_scenario] * 2,
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
    if indicator == 'Secchi':
        # y_range = [0, y_max_plot]

        # if y_median_scenario:
        #     if y_median_scenario > ref_median:
        #         y_range = [0, y_median_scenario * 1.1]
        units = 'm'
    else:
        units = 'mg/m³'

    y_range = [0, y_max]

    # print(y_range)

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
    fig.update_yaxes(range=y_range)

    if measured:
        fig.update_xaxes(range=[start_x, end_x])
    else:
        fig.update_xaxes(range=[start_x, end_x], visible=False)

    return fig


def make_graph(fig=None):
    """
    I need to make a figure this way because of a bug in the reset axes of the figure.
    """
    if fig is not None:
        return html.Div(dcc.Graph(id='ts_plot', figure=fig, config={'responsive': True, 'displaylogo': False, 'modeBarButtonsToRemove': ['zoomIn', 'zoomOut', 'select2d', 'lasso2d']}, style={'height': '50vh'}), id='ts_plot_div', style={'width': '100%', 'height': 'auto', 'margin': "auto", "display": "block"})
    else:
        return html.Div(id='ts_plot_div', style={'width': '100%', 'height': 'auto', 'margin': "auto", "display": "block"})








































































