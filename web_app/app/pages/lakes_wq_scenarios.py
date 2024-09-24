#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 13:37:46 2022

@author: mike
"""
import os
import dash
from dash import dcc, html, dash_table, callback, ctx
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import dash_leaflet as dl
import dash_leaflet.express as dlx
from dash_extensions.javascript import assign, arrow_function
import pandas as pd
import numpy as np
from copy import deepcopy
import base64
import geobuf
import booklet
import npsfm
import plotly.graph_objects as go
# import plotly.express as px
from dash_iconify import DashIconify
import contextily as cx
from pylatex import Document, Section, Subsection, Command, Figure, LongTable, MultiColumn, Tabular, MultiRow, LongTabu, Tabularx, Package, ColumnType, Table
from pylatex.utils import italic, NoEscape
import matplotlib.pyplot as plt
import matplotlib
import tempfile
import geopandas as gpd
import orjson
from shapely import Point
from scipy.stats import norm

# from .app import app
# from . import utils

# from app import app
# import utils

from utils import parameters as param
# from utils import components as gc
from utils import utils

# Set Agg on matplotlib for pylatex to work
matplotlib.use('Agg')

##########################################
### Parameters


dash.register_page(
    __name__,
    path='/lakes-wq-scenarios',
    title='Water Quality Scenarios',
    name='lakes_wq_scenarios',
    description='Lakes Water Quality Scenarios'
)

### Handles
catch_style_handle = assign("""function style(feature) {
    return {
        fillColor: 'grey',
        weight: 2,
        opacity: 1,
        color: 'black',
        fillOpacity: 0.1
    };
}""", name='rivers_catch_style_handle_sites')

base_reach_style_handle = assign("""function style3(feature) {
    return {
        weight: 2,
        opacity: 0.75,
        color: 'grey',
    };
}""", name='rivers_base_reach_style_handle_sites')

sites_points_handle = assign("""function rivers_sites_points_handle(feature, latlng, context){
    const {classes, colorscale, circleOptions, colorProp} = context.props.hideout;  // get props from hideout
    const value = feature.properties[colorProp];  // get value the determines the fillColor
    for (let i = 0; i < classes.length; ++i) {
        if (value == classes[i]) {
            circleOptions.fillColor = colorscale[i];  // set the color according to the class
        }
    }

    return L.circleMarker(latlng, circleOptions);
}""", name='rivers_sites_points_handle_sites')

draw_marae = assign("""function(feature, latlng){
const flag = L.icon({iconUrl: '/assets/nzta-marae.svg', iconSize: [20, 30]});
return L.marker(latlng, {icon: flag});
}""", name='rivers_sites_marae_handle')

lake_poly_style_handle = assign("""function style(feature) {
    return {
        fillColor: 'rgb(166,189,219)',
        weight: 2,
        opacity: 1,
        color: 'black',
        fillOpacity: 1
    };
}""", name='lake_poly_style_handle')

# lake_id = 11133
# indicator = 'TN'
# indicator = 'Secchi'
# lake_name = 'Lake Rotorua'
# lake_data = {'name': 'Lake Rotorua', 'residence_time': 886.1516098743485, 'max_depth': 44.79999923706055, 'mean_depth': 14.93333275968287, 'p_residence_time': 0.9674990155143122, 'n_residence_time': 1, 'regional_council': 'Bay of Plenty', 'area_ha': 2000, 'cluster': False}

# lake_id = 54734
# indicator = 'TN'
# indicator = 'Secchi'
# lake_name = 'Lake Taupo (Taupomoana)'
# lake_data = {'name': 'Lake Taupo                    (Taupomoana)', 'residence_time': 3506.1802279694543, 'max_depth': 162.8000030517578, 'mean_depth': 54.266667837818886, 'p_residence_time': 0.9833922927018965, 'n_residence_time': 1, 'regional_council': 'Waikato', 'area_ha': 2000, 'cluster': False}

# lake_id = 47579
# indicator = 'TN'

# lake_id = 50782
# indicator = 'TN'
# lake_name = 'Lake Waikare'
# lake_data = {'name': 'Lake Waikare', 'residence_time': 81.68687349644028, 'max_depth': 1.7999999523162842, 'mean_depth': 0.5999999759085639, 'p_residence_time': 0.900379346306269, 'n_residence_time': 1, 'area_ha': 3437.4248046875, 'regional_council': 'Waikato', 'cluster': False}


###############################################
### Initial processing

# with booklet.open(eco_sites_path, 'r') as f:
#     catches = [int(c) for c in f]



# catches.sort()
# indicators = list(param.rivers_indicator_dict.keys())
# indicators.sort()

###############################################
### App layout


def layout():
    layout = dmc.Container(
        fluid=True,
        # size='xl',
        # px=0,
        # py=0,
        # my=0,
        # mx=0,
        # ml=0,
        # pl=0,
        # style={
        #     'margin': 0,
        #     'padding': 0
        #     },
        children=[
            dmc.Grid(
                columns=7,
                # grow=True,
                # justify="space-between",
                # align='stretch',
                # style={
                #     'margin': 0,
                #     'padding': 0
                #     },
                children=[
                    dmc.Col(
                        span=3,
                        # style={
                        #     # 'margin-top': 20,
                        #     'margin': 0,
                        #     'padding': 0,
                        #     },
                        children=[
                            html.Label('(1) Select a lake on the map:', style={'font-size': param.header_font_size}),
                            dmc.Text(id='lake_name', weight=700, style={'margin-top': 10}),
                            # dmc.Text('(2) Select an indicator:', style={'margin-top': param.space_between_content, 'font-size': param.header_font_size}),
                            # dcc.Dropdown(options=[{'label': val, 'value': key} for key, val in param.lakes_indicator_dict.items()], id='indicator', optionHeight=40, clearable=False, style={'margin-bottom': param.space_between_content}),
                            dmc.Group(
                                [
                                    html.Label('(2) Input catchment scenario data:', style={'font-size': param.header_font_size}
                                        ),
                                    dmc.HoverCard(
                                        position='left',
                                        withArrow=True,
                                        width=250,
                                        shadow="md",
                                        openDelay=1000,
                                        children=[
                                            dmc.HoverCardTarget(DashIconify(icon="material-symbols:help", width=25)),
                                            dmc.HoverCardDropdown(
                                                dcc.Markdown(
                                                    """
                                                    **Press enter after input to confirm**. The land cover **Native Vegetation** cannot be improved/mitigated. The land cover **Other** includes areas like lakes and river beds that cannot be changed in any way.
                                                    """,
                                                    style={
                                                        'font-size': 14,
                                                        }
                                                    # size="sm",
                                                )
                                            ),
                                        ],
                                    ),
                                    ],
                                style={'margin-top': param.space_between_content}
                                ),
                            dash_table.DataTable(
                                data=[],
                                style_table={'overflowX': 'auto'},
                                style_cell={
                                    'whiteSpace': 'normal',
                                    'height': 'auto',
                                    'minWidth': '50px',
                                    'maxWidth': '100px',
                                    'width': '80px',
                                    'font-size': 12,
                                    'font-family':'sans-serif'
                                },
                                columns=[{'name': name, 'id': key, 'editable': (key in ('phosphorus', 'nitrogen', 'new_land_area'))} for key, name in param.tbl_cols_dict.items()],
                                id='lc_tbl',
                                merge_duplicate_headers=True,
                                style_header={
                                    'textAlign': 'center'},
                                style_cell_conditional=[
                                    {'if': {'column_id': 'land_cover'},
                                      'width': '90px',
                                      'textAlign': 'left'},
                                    {'if': {'column_id': 'land_area'},
                                      'width': '50px'},
                                    {'if': {'column_id': 'new_land_area'},
                                      'width': '50px'},
                                ],
                                style_header_conditional=[{
                                    'if': {'column_id': 'nitrogen'},
                                    'font-weight': 'bold'
                                    },
                                    {
                                        'if': {'column_id': 'phosphorus'},
                                        'font-weight': 'bold'
                                    },
                                    {
                                    'if': {
                                        'column_id': 'new_land_area'},
                                    'font-weight': 'bold'
                                }
                                ]),
                                dmc.Button('Reset to defaults', id='tbl_reset_btn', color='red', variant="light"),
                                dmc.Button('Run catchment scenario', id='tbl_calc_btn', color='green', variant="light"),
                                dcc.Markdown('', style={
                                    'margin-top': 10,
                                    'margin-bottom': param.space_between_content,
                                    'textAlign': 'left',
                                                }, id='tbl_error_text'),

                                dmc.Group(
                                    [
                                        html.Label('(3) Select indicator for results:', style={'font-size': param.header_font_size}
                                            ),
                                        dmc.HoverCard(
                                            position='bottom',
                                            withArrow=True,
                                            width=250,
                                            shadow="md",
                                            openDelay=1000,
                                            children=[
                                                dmc.HoverCardTarget(DashIconify(icon="material-symbols:help", width=25)),
                                                dmc.HoverCardDropdown(
                                                    dcc.Markdown(
                                                        """
                                                        **TN**: Total Nitrogen, **TP**: Total Phosphorus, **Chla**: Chlorophyll a, **Secchi**: Secchi Depth
                                                        """,
                                                        style={
                                                            'font-size': 14,
                                                            }
                                                        # size="sm",
                                                    )
                                                ),
                                            ],
                                        ),
                                        ],
                                    # style={'margin-top': param.space_between_content}
                                    ),
                                dmc.Tabs(
                                    [
                                        dmc.TabsList(
                                            [
                                                dmc.Tab('TN', value='TN'),
                                                dmc.Tab('TP', value='TP'),
                                                dmc.Tab('Chla', value='CHLA'),
                                                dmc.Tab('Secchi', value='Secchi'),
                                            ]
                                        ),
                                    ],
                                    id="results_tabs",
                                    value='TN',
                                ),
                                html.Div(id="stats_tbl",
                                         # style={"paddingTop": 10}
                                         ),
                                dcc.Loading(
                                type="default",
                                children=[
                                    dmc.Button('Download pdf report', color='blue', variant="light", id='dl_btn', disabled=True, style={'margin-top': 10}),
                                    dcc.Download(id="dl_pdf"),
                                    ],
                                ),
                        ]
                        ),
                    dmc.Col(
                        span=4,
                        # style={
                        #     'margin-top': 20
                        #     },
                        children=dmc.Stack(
                            [html.Div([
                            dl.Map(center=param.center, zoom=param.zoom, children=[
                                dl.LayersControl([
                                    dl.BaseLayer(dl.TileLayer(attribution=param.attribution, opacity=0.7), checked=True, name='OpenStreetMap'),
                                    dl.BaseLayer(dl.TileLayer(url='https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png', attribution='Map data: © OpenStreetMap contributors, SRTM | Map style: © OpenTopoMap (CC-BY-SA)', opacity=0.6), checked=False, name='OpenTopoMap'),
                                    dl.BaseLayer(dl.TileLayer(url='https://basemaps.linz.govt.nz/v1/tiles/aerial/WebMercatorQuad/{z}/{x}/{y}.webp?api=d01hgvvhtwhnhestm1bdf6mtvr5', attribution='Map data: © <a href="//www.linz.govt.nz/linz-copyright">LINZ CC BY 4.0</a> © <a href="//www.linz.govt.nz/data/linz-data/linz-basemaps/data-attribution">Imagery Basemap contributors</a>', opacity=0.6), checked=False, name='Aerial Imagery'),
                                    dl.Overlay(dl.LayerGroup(dl.GeoJSON(data='', format="geobuf", id='catch_map', zoomToBoundsOnClick=True, zoomToBounds=True, options=dict(style=catch_style_handle))), name='Catchments', checked=True),
                                    # dl.Overlay(dl.LayerGroup(dl.GeoJSON(data='', format="geobuf", id='marae_map_sites', zoomToBoundsOnClick=False, zoomToBounds=False, options=dict(pointToLayer=draw_marae))), name='Marae', checked=False),
                                    dl.Overlay(dl.LayerGroup(dl.GeoJSON(data='', format="geobuf", id='reach_map', options=dict(style=base_reach_style_handle), hideout={})), name='Rivers', checked=True),
                                    dl.Overlay(dl.LayerGroup(dl.GeoJSON(data='', format="geobuf", id='lake_poly_map', zoomToBoundsOnClick=True, cluster=False, options=dict(style=lake_poly_style_handle))), name='Lakes', checked=True),
                                    dl.Overlay(dl.LayerGroup(dl.GeoJSON(url=str(param.lakes_points_path), format="geobuf", id='lakes_points_map', zoomToBoundsOnClick=True, cluster=True)), name='Lakes Points', checked=True),
                                    ],
                                    id='layers'
                                    ),
                                # gc.colorbar_power,
                                # html.Div(id='colorbar', children=colorbar_base),
                                # dmc.Group(id='colorbar', children=colorbar_base),
                                # dcc.Markdown(id="info_sites", className="info", style={"position": "absolute", "top": "10px", "right": "160px", "z-index": "1000"})
                                ],
                                style={'width': '100%', 'height': param.map_height, 'margin': "auto", "display": "block"}
                                ),

                            ],
                            # className='five columns', style={'margin': 10}
                            ),
                            # html.Div(id='ts_plot',
                            #          style={'width': '100%', 'height': '40vh', 'margin': "auto", "display": "block"}),
                            utils.make_graph()
                            ]
                            ),
                        ),
                    ]
                    ),
            dcc.Store(id='lake_id', data=0),
            dcc.Store(id='lake_data', data=None),
            dcc.Store(id='calc_ready', data=0),
            dcc.Store(id='stats', data={}),
            dcc.Store(id='conc_factors', data=None),
            dcc.Store(id='box_plot_fig', data=None),
            ]
        )

    return layout


###############################################
### Callbacks


@callback(
    Output('lake_id', 'data'),
    Output('lake_data', 'data'),
    Output('lake_name', 'children'),
    [Input('lakes_points_map', 'click_feature')]
    )
def update_lake_id(feature):
    """

    """
    # print(ds_id)
    lake_id = 0
    lake_name = ''
    lake_data = None
    if feature is not None:
        # print(feature)
        if not feature['properties']['cluster']:
            lake_id = int(feature['id'])
            lake_data = feature['properties']
            lake_name = lake_data['name']
            # print(lake_id)
            # print(feature)

    return lake_id, lake_data, lake_name


@callback(
        Output('lake_poly_map', 'data'),
        Input('reach_map', 'data'),
        State('lake_id', 'data'),
        )
# @cache.memoize()
def update_lake(reach_map, lake_id):
    if lake_id > 0:
        # with booklet.open(param.lakes_poly_gbuf_path, 'r') as f:
        #     data = base64.b64encode(f[int(lake_id)]).decode()

        data = base64.b64encode(utils.get_value(param.lakes_poly_path, lake_id)).decode()
    else:
        data = ''

    return data


@callback(
        Output('catch_map', 'data'),
        Input('lake_id', 'data'),
        prevent_initial_call=True
        )
# @cache.memoize()
def update_catchment(lake_id):
    if lake_id > 0:
        data = base64.b64encode(utils.get_value(param.lakes_catch_path, lake_id)).decode()

    else:
        data = ''

    return data


@callback(
        Output('reach_map', 'data'),
        Input('lake_id', 'data'),
        prevent_initial_call=True
        )
# @cache.memoize()
def update_reaches(lake_id):
    if lake_id > 0:
        data = base64.b64encode(utils.get_value(param.lakes_catch_reaches_path, lake_id)).decode()

    else:
        data = ''

    return data


@callback(
        Output('calc_ready', 'data'),
        Output('tbl_error_text', 'children'),
        Input('lc_tbl', 'data'),
        prevent_initial_call=True
        )
def check_tbl_values(lc_tbl):
    """

    """
    if lc_tbl:
        try:
            data = pd.DataFrame(lc_tbl).set_index('land_cover').astype('int16')
        except ValueError:
            return 0, 'All inputs must be integers'
        # print(data)

        lcs = data.index

        mitigation_cols = ['phosphorus', 'nitrogen']

        if 'Other' in lcs:
            unchange = data.loc['Other']
            if unchange['land_area'] != unchange['new_land_area']:
                return 0, '**Other** land areas cannot be changed'
            if (unchange[mitigation_cols] != 0).any():
                return 0, '**Other** lands cannot be mitigated'

        if 'Native Vegetation' in lcs:
            unimprove = data.loc['Native Vegetation']
            if (unimprove[mitigation_cols] != 0).any():
                return 0, '**Native Vegetation** lands cannot be mitigated'

        if (data[mitigation_cols] < 0).any(axis=1).any():
            return 0, "Land mitigation %'s must be >= 0 and <= 100"
        if (data[mitigation_cols] > 100).any(axis=1).any():
            return 0, "Land mitigation %'s must be >= 0 and <= 100"
        if any(data['new_land_area'] < 0):
            return 0, "New land area %'s must be >= 0"
        if data['new_land_area'].sum() != 100:
            leftover = 100 - data['new_land_area'].sum()
            if leftover > 0:
                text = f"New land area % must add up to 100%, please allocate {leftover}%"
            else:
                text = f"New land area % must add up to 100%, please remove {-leftover}%"
            return 0, text

        return 1, ''
    else:
        return 0, ''


@callback(
        Output('lc_tbl', 'data'),
        Output('conc_factors', 'data'),
        Input('lake_id', 'data'),
        # Input('indicator', 'value'),
        Input('tbl_reset_btn', 'n_clicks'),
        Input('tbl_calc_btn', 'n_clicks'),
        # State('nzsegment', 'data'),
        State('lc_tbl', 'data'),
        State('calc_ready', 'data'),
        # State('lake_data', 'data'),
        prevent_initial_call=True
        )
def update_mitigation_tbl(lake_id, reset, calc, lc_tbl, calc_ready):
    """

    """
    trig = ctx.triggered_id

    # print(indicator)

    if lake_id > 0:
        lc_tbl, conc_factors = utils.calc_scenario_tbl(lake_id, lc_tbl, calc_ready, trig)

    else:
        conc_factors = None
        lc_tbl = []

    return lc_tbl, conc_factors


@callback(
        Output('stats', 'data'),
        Output('dl_btn', 'disabled'),
        Input('conc_factors', 'data'),
        State('lake_id', 'data'),
        # State('indicator', 'value'),
        State('stats', 'data'),
        State('lake_data', 'data'),
        prevent_initial_call=True
        )
def calc_stats(conc_factors, lake_id, stats, lake_data):
    # trig = ctx.triggered_id

    # print(stats_tbl)
    # print(source)

    disabled = True

    if lake_id > 0:

        ### nps-fm stats
        stats = {}

        for indicator in param.lakes_indicator_dict:
            ind_stats = {}
            if indicator in param.nps_mapping:
                results_str1 = param.indicator_str_format[indicator]
                nps = npsfm.NPSFM(param.assets_path)
                nps_param = param.nps_mapping[indicator]
    
                limits = nps.add_limits('lake', nps_param, lfenzid=lake_id)
                # len_limits = len(limits)
                initial_conc = 0
                for i, band in enumerate(limits):
                    median = limits[band]['median'][1]
                    if median >= 10000:
                        band_str = f'> {initial_conc}'
                        ind_stats[f'Band {band}'] = band_str
                        break
                    else:
                        band_str = f'{initial_conc} - {median}'
                        ind_stats[f'Band {band}'] = band_str
    
                    initial_conc = median
    
                bl_limit = nps.bottom_line_limit['median'][1]
                if bl_limit >= 10000:
                    bl_str = 'No bottom line'
                else:
                    bl_str = results_str1.format(bl_limit)
    
                ind_stats['Bottom line'] = bl_str

            stats[indicator] = ind_stats

        ## Ref conc
        # ref_conc0 = utils.get_value(param.lakes_ref_conc_path, lake_id)

        # for indicator in param.lakes_indicator_dict:
        #     results_str1 = param.indicator_str_format[indicator]
        #     ref_conc = ref_conc0[indicator]
        #     stats[indicator].append({'name': 'Reference', 'conc': results_str1.format(ref_conc)})

        ## Current - measured
        try:
            stats0 = utils.get_value(param.lakes_moni_medians_path, lake_id)

            for indicator in param.lakes_indicator_dict:
                if indicator not in stats0:
                    raise ValueError()

            for indicator in param.lakes_indicator_dict:
                results_str1 = param.indicator_str_format[indicator]
                meas_median = stats0[indicator]

                stats[indicator]['Current (measured)'] = results_str1.format(meas_median)

        except:
            ## Current - modelled
            stats0 = utils.get_value(param.lakes_model_medians_path, lake_id)
            for indicator in param.lakes_indicator_dict:
                results_str1 = param.indicator_str_format[indicator]
                model_median = stats0[indicator]
        
                stats[indicator]['Current (modelled)'] = results_str1.format(model_median)

        # ## Achievable concs
        # stats0 = utils.get_value(param.lake_achievable_conc_path, lake_id)


        ### Scenario calcs
        if conc_factors is not None:
        # if trig == 'tbl_calc_btn':
            region = lake_data['regional_council']
            t = lake_data['residence_time']
            z_max = lake_data['max_depth']
            scenario_stats = utils.est_ind_scenario_conc(stats0, conc_factors, region, t, z_max)
            # print(conc_factors)
            # print(stats0)
            # print(scenario_stats)
            for indicator in param.lakes_indicator_dict:
                results_str1 = param.indicator_str_format[indicator]
                scenario_stat = scenario_stats[indicator]
                stats[indicator]['Scenario'] = results_str1.format(scenario_stat)

            disabled = False
        else:
            for indicator in param.lakes_indicator_dict:
                stats[indicator]['Scenario'] = 'Press the Run catchment scenario button'

    else:
        stats = {}

    # print(stats)
    # print(nzsegment)

    return stats, disabled


@callback(
        Output('stats_tbl', 'children'),
        Input('results_tabs', 'value'),
        Input('stats', 'data'),
        # State('indicator', 'value'),
        prevent_initial_call=True
        )
def update_stats_table(tab, stats):

    stats_tbl_data = []
    # print(tab)
    if tab in stats:
        data = stats[tab]
        for k, v in data.items():
            if ('Reference' in k) or ('Current' in k) or ('Scenario' in k):
                stats_tbl_data.append({'name': k, 'conc': v})

    stats_tbl = utils.make_results_table(stats_tbl_data, tab)

    return stats_tbl


@callback(
        Output('ts_plot_div', 'children'),
        Input('results_tabs', 'value'),
        Input('stats', 'data'),
        State('lake_id', 'data'),
        prevent_initial_call=True
        )
def update_box_plot(tab, stats, lake_id):

    # print(stats_tbl)
    if tab in stats:
        data_tbl = stats[tab]
        fig = utils.make_fig_ind(data_tbl, lake_id, tab)
    else:
        fig = None

    return utils.make_graph(fig)


@callback(
        Output('dl_pdf', 'data'),
        Input("dl_btn", "n_clicks"),
        # State('ts_plot', 'figure'),
        State('lc_tbl', 'data'),
        State('stats', 'data'),
        State('lake_id', 'data'),
        State('lake_name', 'children'),
        # State('nzsegment', 'data'),
        State('conc_factors', 'data'),
        State('lake_data', 'data'),
        prevent_initial_call=True,
        )
def make_pdf_report(n_clicks, lc_tbl, stats, lake_id, lake_name, conc_factors, lake_data):
    """

    """
    ## Pre-calcs
    # Determine whether measured or modelled
    measured = True
    measured_text = ' (measured)'
    for k in stats['TN']:
        if 'modelled' in k:
            measured = False
            measured_text = ' (modelled)'

    current_text = 'Current' + measured_text
    scenario_text = 'Scenario'

    current_scenario = {}
    improve_perc_dict = {}
    improve_text_dict = {}
    for ind, res in stats.items():
        ind_stats = {k.split(' (')[0]: float(v) for k, v in res.items() if k in (current_text, scenario_text, 'Reference')}
        current_scenario[ind] = ind_stats

        improve_perc = int(round((1 - (ind_stats['Scenario']/ind_stats['Current'])) * 100))

        if ind == 'Secchi':
            if improve_perc < 0:
                improve_perc = improve_perc * -1
                improve_text = 'improvement'
            else:
                improve_text = 'degredation'
        else:
            if improve_perc < 0:
                improve_perc = improve_perc * -1
                improve_text = 'degredation'
            else:
                improve_text = 'improvement'

        improve_perc_dict[ind] = improve_perc
        improve_text_dict[ind] = improve_text

    ## Reference conc
    ref_conc0 = utils.get_value(param.lakes_ref_conc_path, lake_id)

    ## Conc tbl
    row_names = ['Band A', 'Band B', 'Band C', 'Band D', 'Bottom line']
    if measured:
        row_names += ['Current (measured)', 'Scenario']
    else:
        row_names += ['Current (modelled)', 'Scenario']

    conc_tbl_data = []
    for row_name in row_names:
        row = [row_name]
        for ind in param.lakes_indicator_dict:
            val = None
            for k, data in stats[ind].items():
                if k == row_name:
                    row.append(data)
                    val = data
                    break
            if val is None:
                row.append('')
        conc_tbl_data.append(row)

    ## Measured data
    if measured:
        stats0 = utils.get_value(param.lakes_moni_data_path, lake_id)
    
        grp1 = stats0.reset_index().groupby('indicator')
    
        # obs_count = grp1['value'].count().to_dict()
        obs_min_dates = grp1['date'].first().to_dict()
        obs_max_dates = grp1['date'].last().to_dict()
        obs_medians = grp1['value'].median().to_dict()

    # scenario_conc = current_scenario['Scenario']
    # ref_mean, ref_lower, ref_upper = utils.sep_reference_values([s['conc'] for s in stats if s['name'] == 'Reference'][0])
    # nps_check = [True for s in stats if s['name'] == 'Bottom line']

    lc_data = utils.get_value(param.lakes_lc_red_yields_path, lake_id).sort_index()
    lc_data1 = lc_data[[col for col in lc_data.columns if 'reduction' not in col]].round(2).reset_index().copy()
    lc_data1['area_ha'] = lc_data1['area_ha'].round().astype(int)
    tot_area = lc_data1['area_ha']['area_ha'].sum()

    # if indicator != 'ECOLI':
    #     ind_name = param.rivers_indicator_dict[indicator].lower()
    # else:
    #     ind_name = param.lakes_indicator_dict[indicator]
    # river_name = site_data['properties']['Catchment']
    agency = lake_data['regional_council']
    lake_area = '{:,.0f}'.format(lake_data['area_ha'])
    residence_time = '{:,.0f}'.format(lake_data['residence_time'])
    max_depth = '{:,.1f}'.format(lake_data['max_depth'])

    # box_plot_fig = utils.decode_obj(box_plot_fig_enc)
    # print(geometry)
    # geom = orjson.loads(geometry.encode())

    ## plot the catchment polygon and river network
    fig, ax = plt.subplots()

    ## Rivers
    lakes_reaches_gbuf = utils.get_value(param.lakes_catch_reaches_path, lake_id)

    lakes_reaches_dict = geobuf.decode(lakes_reaches_gbuf)
    if lakes_reaches_dict is not None:
        lakes_geo = gpd.GeoDataFrame.from_features(lakes_reaches_dict['features'], crs=4326)
        lakes_geo.plot(ax=ax, alpha=0.3, zorder=2)

    ## Catchment
    catch_gbuf = utils.get_value(param.lakes_catch_path, lake_id)

    catch_dict = geobuf.decode(catch_gbuf)
    catch_geo = gpd.GeoDataFrame.from_features(catch_dict['features'], crs=4326)
    catch_geo.plot(ax=ax, color='lightgrey', edgecolor='black', alpha=0.6, zorder=1)

    ## lake polygon location
    poly_gbuf = utils.get_value(param.lakes_poly_path, lake_id)
    poly_dict = geobuf.decode(poly_gbuf)
    poly_geo = gpd.GeoDataFrame.from_features(poly_dict['features'], crs=4326)
    poly_geo.plot(ax=ax, color='lightblue', edgecolor='black', alpha=1, zorder=3)

    # point = gpd.GeoSeries([Point(lake_data['geometry']['coordinates'])], crs=4326)
    # point.plot(ax=ax, color='black')
    # ax.set_xlim(x_min - b, x_max + b)
    # ax.set_ylim(y_min - b, y_max + b)
    # ax.set_extent((x_min - b, y_min - b, x_max + b, y_max + b))
    cx.add_basemap(ax, crs=catch_geo.crs, source=cx.providers.OpenStreetMap['Mapnik'])

    fig.tight_layout()

    ## Make the document
    lc_tbl_dict = {k: l[1] for k, l in param.tbl_cols_dict.items()}

    with tempfile.TemporaryDirectory() as path:
        doc = Document(os.path.join(path, f'lake_{lake_id}'), geometry_options={"right": "2cm", "left": "2cm"})
        doc.packages.append(Package('array'))
        doc.packages.append(Package('datetime2'))
        doc.packages.append(Package('float'))
        # right_fixed_col = ColumnType('R', 'p{#1}', '\raggedright')
        doc.preamble.append(NoEscape(r'\newcolumntype{R}[1]{>{\raggedleft\arraybackslash}p{#1}}'))
        doc.preamble.append(NoEscape(r'\newcolumntype{L}[1]{>{\raggedright\arraybackslash}p{#1}}'))
        doc.preamble.append(NoEscape(r'\newcolumntype{C}[1]{>{\centering\arraybackslash}p{#1}}'))

        if lake_name == 'No name':
            doc.preamble.append(Command('title', f'Scenario builder results for the lake with the DOC lake id {lake_id}'))
        else:
            doc.preamble.append(Command('title', f'Scenario builder results for {lake_name}'))
        # doc.preamble.append(catch_name)
        # doc.preamble.append(Command('date', NoEscape(r'\today')))
        # doc.preamble.append(NoEscape(r'\today'))
        doc.append(NoEscape(r'\maketitle'))

        ## Catchment location section
        with doc.create(Section('Monitoring site and catchment information')):
            if lake_name == 'No name':
                text = f'The user-selected lake has no official name and is within {agency}. '
            else:
                text = f'The user-selected lake is named "{lake_name}" and is within {agency}. '

            text += f'The lake is also given the ID of {lake_id} by DOC. The lake has an estimated area of {lake_area} ha, a maximum depth of {max_depth} m, and a water residence time of {residence_time} years. The surrounding surface water catchment has an estimated area of {tot_area:,} ha '

            # if catch_name != 'Unnamed':
            #     text += f'It is located within the larger {catch_name} catchment '

            doc.append(text)
            doc.append(NoEscape(r'(Figure \ref{fig:catchment}).'))
            # position='htbp'
            with doc.create(Figure(position='H')) as plot:
                plot.add_plot(width=NoEscape(r'0.8\textwidth'))
                plot.add_caption('The lake and the associated surface water catchment. ')
                plot.append(Command('label', 'fig:catchment'))

        ## Section 1
        sec1 =  Section('Land area contributions and mitigations')
        doc.append(sec1)

        sec_text = r"The current and scenario land areas and mitigations are shown in Table \ref{tab:mitigations}. "
        sec_text += "The nitrogen and phosphorus mitigations are used to estimate the improvements across all indicators. "

        doc.append(NoEscape(sec_text))

        ## Land use table
        base_tbl = Table(position="htb")
        base_tbl.add_caption('Land area contributions and mitigations')
        base_tbl.append(NoEscape(r'\centering'))
        data_table = Tabular('| l | R{0.06\linewidth} | R{0.08\linewidth} | R{0.10\linewidth} | R{0.12\linewidth} | R{0.12\linewidth} | R{0.06\linewidth} |')
        data_table.add_hline()
        data_table.add_row(['', MultiColumn(3, align='c|', data='Current'), MultiColumn(3, align='c|', data='Scenario')], strict=False)
        data_table.add_row(list(lc_tbl_dict.values()), strict=False)
        data_table.add_hline()
        for values in lc_tbl:
            row = [str(int(val)) if isinstance(val, (float, int)) else val for val in values.values()]
            data_table.add_row(row, strict=False)
        data_table.add_hline()
        base_tbl.append(data_table)
        # base_tbl.add_caption('Land areas and mitigations')
        base_tbl.append(Command('label', 'tab:mitigations'))

        doc.append(base_tbl)

        table_mitigations = r'\ref{tab:mitigations}'

        ## Section 2
        # sec2 =  Section('Land areas and mitigations')
        sec2 =  Section('Current estimates')
        doc.append(sec2)

        header = ['Name', 'Total nitrogen (mg/m³)', 'Total phosphorus (mg/m³)', 'Chlorophyll a (mg/m³)', 'Secchi depth (m)']
        base_tbl2 = Table(position="htb")
        base_tbl2.add_caption('Site median concentrations')
        base_tbl2.append(NoEscape(r'\centering'))
        data_table2 = Tabular("| l | R{0.09\linewidth} | R{0.11\linewidth} | R{0.12\linewidth} | R{0.06\linewidth} |")
        data_table2.add_hline()
        data_table2.add_row(header, strict=False)
        data_table2.add_hline()
        for row in conc_tbl_data:
            # row = [values['name'], values['conc']]
            data_table2.add_row(row, strict=False)
        data_table2.add_hline()
        base_tbl2.append(data_table2)
        # base_tbl2.add_caption('Site median concentrations')
        base_tbl2.append(Command('label', 'tab:conc'))

        doc.append(base_tbl2)

        table_conc_ref = r'\ref{tab:conc}'

        if measured:
            sec_text = f'This lake has monitoring data and subsequently the median concentration estimates are from the measured data (Table {table_conc_ref}). '
        else:
            sec_text = f'This lake does not have monitoring data and subsequently the median concentration estimates are from modelled results (Table {table_conc_ref}). '

        doc.append(NoEscape(sec_text))

        # 2.1 - Indicators
        # sec2a = Subsection('Indicators')
        for ind, res in current_scenario.items():
            # rounding = param.indicator_rounding[ind]
            rounding_str = param.indicator_str_format[ind]
            if ind == 'Secchi':
                units = 'm'
            else:
                units = 'mg/m³'

            ind_name = param.lakes_indicator_dict[ind]
            # sec2a = Subsection(f'{ind_name}')

            if measured:
                obs_min_date = str(obs_min_dates[ind])
                obs_max_date = str(obs_max_dates[ind])
                obs_median = rounding_str.format(obs_medians[ind])

                sec_text = f'The {ind_name.lower()} current measured estimate of the median concentration for the 5 year period from an initial observation at {obs_min_date} to a final observation at {obs_max_date} is {obs_median} {units}. '
            else:
                model_median = rounding_str.format(res['Current'])
                sec_text = f'The {ind_name.lower()} current modelled estimate of the median concentration is {model_median} {units}. '

            # table_conc_ref = r'\ref{tab:conc}'

            # sec_text += f'(Table {table_conc_ref}). '

            ref_median = rounding_str.format(ref_conc0[ind])

            sec_text += f'The estimated reference median concentration is {ref_median} {units}. '

            doc.append(NoEscape(sec_text))

        sec_text = f'The relevant NPS-FM 2020 bands associated with the median concentration attribute states are also listed in Table {table_conc_ref}. The reference concentration is an independent estimate of what the natural concentration of the lake would have been before human impact. '

        doc.append(NoEscape(sec_text))

        ## Section 3
        sec3 = Section('Scenario outcomes')
        doc.append(sec3)

        sec_text = ''

        for ind, res in current_scenario.items():
            rounding_str = param.indicator_str_format[ind]

            if ind == 'Secchi':
                units = 'm'
            else:
                units = 'mg/m³'

            ind_name = param.lakes_indicator_dict[ind]
            # sec3a = Subsection(f'{ind_name}')

            # table_conc_ref = r'\ref{tab:conc}'

            improve_text = improve_text_dict[ind]
            improve_perc = improve_perc_dict[ind]
            scenario_conc = rounding_str.format(res['Scenario'])

            sec_text += f'The {ind_name.lower()} estimated {improve_text} at the lake is {improve_perc}\% with a resulting median concentration of {scenario_conc} {units}. '

            # sec_text += f'the estimated {improve_text} at the lake is {improve_perc}\% with a resulting median concentration of {scenario_conc} {units} '
            # doc.append(NoEscape(sec_text))

        sec_text += f'All of the scenario results were estimated using the user-defined mitigations/land use changes as shown in Table {table_mitigations}. '
        sec_text += 'Refer to the user guide for more details on the background and methodology of the calculations.'

        # sec3a.append(NoEscape(sec_text))
        doc.append(NoEscape(sec_text))

        doc.generate_pdf(clean_tex=True)

        pdf_path = os.path.join(path, f'lake_{lake_id}.pdf')
        with open(pdf_path, 'rb') as f:
            # pdf_bytes = io.BytesIO(f.read())
            pdf_bytes = f.read()

        os.remove(pdf_path)

    return dcc.send_bytes(pdf_bytes, f'lake_{lake_id}.pdf', 'pdf')








