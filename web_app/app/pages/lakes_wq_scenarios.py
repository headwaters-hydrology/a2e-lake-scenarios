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
    path='/rivers-wq-scenarios',
    title='Water Quality Scenarios',
    name='rivers_wq_scenarios',
    description='River Water Quality Scenarios'
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



# site_id = 'ebop-00010'
# indicator = 'DRP'
# indicator = 'ECOLI'
# nzsegment = '4084132'
# site_name = 'Ngongotaha at SH36'
# source = 'CD/L'
# site_data = {'type': 'Feature', 'geometry': {'type': 'Point', 'coordinates': [174.549643, -39.754258]}, 'id': 'trc-00050', 'properties': {'LawaSiteID': 'trc-00050', 'Region': 'taranaki', 'Agency': 'Taranaki Regional Council', 'SiteID': 'Whenuakura at Nicholson Road', 'CouncilSiteID': 'wnr000450', 'Catchment': 'Whenuakura', 'Landuse': 'Rural', 'Altitude': 'Lowland', 'nzsegment': 6213881, 'tooltip': 'Whenuakura at Nicholson Road', 'cluster': False}}

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
                            dmc.Text('(2) Select an indicator:', style={'margin-top': param.space_between_content, 'font-size': param.header_font_size}),
                            dcc.Dropdown(options=[], id='indicator', optionHeight=40, clearable=False, style={'margin-bottom': param.space_between_content}),
                            dmc.Group(
                                [
                                    html.Label('(3) Input catchment scenario data:', style={'font-size': param.header_font_size}
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
                                columns=[{'name': name, 'id': key, 'editable': (key in ('mitigation', 'new_land_area'))} for key, name in param.tbl_cols_dict.items()],
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
                                    'if': {'column_id': 'mitigation_n'},
                                    'font-weight': 'bold'
                                    },
                                    {
                                        'if': {'column_id': 'mitigation_p'},
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
                                        html.Label('(4) Results:', style={'font-size': param.header_font_size}
                                            ),
                                        # dmc.HoverCard(
                                        #     position='bottom',
                                        #     withArrow=True,
                                        #     width=250,
                                        #     shadow="md",
                                        #     openDelay=1000,
                                        #     children=[
                                        #         dmc.HoverCardTarget(DashIconify(icon="material-symbols:help", width=25)),
                                        #         dmc.HoverCardDropdown(
                                        #             dcc.Markdown(
                                        #                 """
                                        #                 **Indicator** shows the results for the selected indicator, while **Periphyton** shows the results of the Q92 Chla estimate based on the selected indicator. Periphyton cannot be estimated from *E. coli*. See the User Guide for shade definition. The Reference values refer to the median, 5th, and 95th percentiles.
                                        #                 """,
                                        #                 style={
                                        #                     'font-size': 14,
                                        #                     }
                                        #                 # size="sm",
                                        #             )
                                        #         ),
                                        #     ],
                                        # ),
                                        ],
                                    # style={'margin-top': param.space_between_content}
                                    ),
                                # dmc.Tabs(
                                #     [
                                #         dmc.TabsList(
                                #             [
                                #                 dmc.Tab('Measured indicator', value='ind'),
                                #                 dmc.Tab('Periphyton', value='peri', disabled=True, id='peri_tab'),
                                #             ]
                                #         ),
                                #     ],
                                #     id="results_tabs",
                                #     value='ind',
                                # ),
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
                                    dl.Overlay(dl.LayerGroup(dl.GeoJSON(data='', format="geobuf", id='catch_map', zoomToBoundsOnClick=False, zoomToBounds=True, options=dict(style=catch_style_handle))), name='Catchments', checked=True),
                                    # dl.Overlay(dl.LayerGroup(dl.GeoJSON(data='', format="geobuf", id='marae_map_sites', zoomToBoundsOnClick=False, zoomToBounds=False, options=dict(pointToLayer=draw_marae))), name='Marae', checked=False),
                                    dl.Overlay(dl.LayerGroup(dl.GeoJSON(data='', format="geobuf", id='reach_map', options=dict(style=base_reach_style_handle), hideout={})), name='Rivers', checked=True),
                                    dl.Overlay(dl.LayerGroup(dl.GeoJSON(data='', format="geobuf", id='lake_poly_map', zoomToBoundsOnClick=True, cluster=False)), name='Lakes', checked=True),
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
                            # dmc.Group(
                            #     children=[
                            #         dmc.Col(html.Div(id='box_plot'), span=6),
                            #         dmc.Col(html.Div(id='bar_plot'), span=6)
                            #         ]
                            #     )

                            # dmc.Group([
                            # dmc.Col(span=2,
                            #         children=html.Div([dcc.Graph(id='box_plot')])
                            #         ),
                            # dmc.Col(span=2,
                            #         children=html.Div([dcc.Graph(id='bar_plot')])
                            #         ),
                            # ]
                            # ),
                            ]
                            ),
                        ),
                    ]
                    ),
            dcc.Store(id='lake_id', data=0),
            # dcc.Store(id='nzsegment', data=''),
            # dcc.Store(id='source', data=''),
            dcc.Store(id='lake_data', data=None),
            dcc.Store(id='calc_ready', data=0),
            dcc.Store(id='stats', data=[]),
            dcc.Store(id='conc_factor', data=1),
            dcc.Store(id='box_plot_fig', data=None),
            ]
        )

    return layout


###############################################
### Callbacks

# @callback(
#     Output('lake_id', 'data'),
#     Output('nzsegment', 'data'),
#     Output('lake_data', 'data'),
#     Output('source', 'data'),
#     [Input('sites_map', 'click_feature')],
#     prevent_initial_call=True
#     )
# def update_lake_id(feature):
#     """

#     """
#     lake_id = ''
#     nzsegment = ''
#     source = ''
#     if feature is not None:
#         if not feature['properties']['cluster']:
#             lake_id = str(feature['id'])
#             nzsegment = str(feature['properties']['nzsegment'])
#             source = feature['properties']['source']
#             # geometry = feature['geometry']

#     # print(feature)
#     # print(source)

#     return lake_id, nzsegment, feature, source


@callback(
    Output('lake_id', 'data'),
    Output('lake_data', 'data'),
    Output('lake_name', 'children')
    [Input('lake_points', 'click_feature')]
    )
def update_lake_id(feature):
    """

    """
    # print(ds_id)
    lake_id = ''
    if feature is not None:
        # print(feature)
        if not feature['properties']['cluster']:
            lake_id = str(feature['id'])
            lake_name = feature['properties']['name']

    return lake_id, feature, lake_name


# @callback(
#     Output('lake_name', 'children'),
#     [Input('lake_id', 'data')],
#     prevent_initial_call=True
#     )
# def update_site_name(lake_id):
#     """

#     """
#     # print(ds_id)
#     if lake_id != '':
#         lake_name = utils.get_value(param.lakes_names_path, lake_id)

#         return lake_name


@callback(
        Output('lake_poly_map', 'data'),
        Input('lake_id', 'data'),
        )
# @cache.memoize()
def update_lake(lake_id):
    if lake_id > 0:
        # with booklet.open(param.lakes_poly_gbuf_path, 'r') as f:
        #     data = base64.b64encode(f[int(lake_id)]).decode()

        data = base64.b64encode(utils.get_value(param.lakes_poly_gbuf_path, lake_id)).decode()
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
        data = base64.b64encode(utils.get_value(param.lakes_catch_gbuf_path, lake_id)).decode()

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
        data = base64.b64encode(utils.get_value(param.lakes_catch_reaches_gbuf_path, lake_id)).decode()

    else:
        data = ''

    return data


@callback(
        Output('indicator', 'options'),
        Output('indicator', 'value'),
        Input('lake_id', 'data'),
        State('indicator', 'value'),
        prevent_initial_call=True
        )
def update_indicators(lake_id, indicator):
    output = []
    if lake_id > 0:
        indicators = utils.get_value(param.lakes_ind_path, lake_id)

        output = [{'label': val, 'value': key} for key, val in param.lakes_indicator_dict.items() if key in indicators]

        if indicator not in indicators:
            indicator = None

    return output, indicator


@callback(
        Output('lc_tbl', 'data'),
        Output('conc_factor', 'data'),
        Input('lake_id', 'data'),
        Input('indicator', 'value'),
        Input('tbl_reset_btn', 'n_clicks'),
        Input('tbl_calc_btn', 'n_clicks'),
        # State('nzsegment', 'data'),
        State('lc_tbl', 'data'),
        State('calc_ready', 'data'),
        prevent_initial_call=True
        )
def update_mitigation_tbl(lake_id, indicator, reset, calc, lc_tbl, calc_ready):
    """

    """
    trig = ctx.triggered_id

    conc_factor = 1

    # print(indicator)

    if (lake_id > 0) and (indicator is not None):
        lc_tbl, conc_factor = utils.calc_scenario_results(lake_id, indicator, lc_tbl, calc_ready, trig)

    else:
        lc_tbl = []

    return lc_tbl, conc_factor


@callback(
        Output('calc_ready', 'data'),
        Output('tbl_error_text', 'children'),
        Input('lc_tbl', 'data'),
        prevent_initial_call=True
        )
def check_tbl_values(tbl_data):
    """

    """
    if tbl_data:
        try:
            data = pd.DataFrame(tbl_data).set_index('land_cover').astype('int16')
        except ValueError:
            return 0, 'All inputs must be integers'
        # print(data)

        lcs = data.index

        if 'Other' in lcs:
            unchange = data.loc['Other']
            if unchange['land_area'] != unchange['new_land_area']:
                return 0, '**Other** land areas cannot be changed'
            if unchange['mitigation'] != 0:
                return 0, '**Other** lands cannot be mitigated'

        if 'Native Vegetation' in lcs:
            unimprove = data.loc['Native Vegetation']
            if unimprove['mitigation'] != 0:
                return 0, '**Native Vegetation** lands cannot be mitigated'

        if any((data['mitigation_n'] < 0) | (data['mitigation_n'] > 100)):
            return 0, "Land mitigation %'s must be >= 0 and <= 100"
        if any((data['mitigation_p'] < 0) | (data['mitigation_p'] > 100)):
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
        Output('stats', 'data'),
        Output('dl_btn', 'disabled'),
        Input('calc_ready', 'data'),
        Input('tbl_calc_btn', 'n_clicks'),
        State('lake_id', 'data'),
        State('indicator', 'value'),
        # State('nzsegment', 'data'),
        # State('source', 'data'),
        State('conc_factor', 'data'),
        State('stats', 'data'),
        prevent_initial_call=True
        )
def calc_stats(calc_ready, calc, lake_id, indicator, conc_factor, stats):
    trig = ctx.triggered_id

    # print(stats_tbl)
    # print(source)

    disabled = True

    if calc_ready == 1:
        if indicator == 'ECOLI':
            # units = 'CFU/100ml'
            results_str1 = '{:.0f}'
            results_str2 = '{:.0f} [ {:.0f} - {:.0f} ]'
            # band_str1 = '{initial_conc} - {median}'
        else:
            # units = 'mg/l'
            results_str1 = '{:.3f}'
            results_str2 = '{:.3f} [ {:.3f} - {:.3f} ]'

        ### Indicator stats
        stats = []

        if indicator in param.nps_mapping:
            nps = npsfm.NPSFM(param.assets_path)
            nps_param = param.nps_mapping[indicator]

            limits = nps.add_limits('lake', nps_param, 1000010)
            # len_limits = len(limits)
            initial_conc = 0
            for i, band in enumerate(limits):
                median = limits[band]['median'][1]
                if median >= 10000:
                    band_str = f'> {initial_conc}'
                    stats.append({'name': f'Band {band}', 'conc': band_str})
                    break
                else:
                    band_str = f'{initial_conc} - {median}'
                    stats.append({'name': f'Band {band}', 'conc': band_str})

                initial_conc = median

            bl_limit = nps.bottom_line_limit['median'][1]
            if bl_limit >= 10000:
                bl_str = 'No bottom line'
            else:
                bl_str = results_str1.format(bl_limit)

            stats.append({'name': 'Bottom line', 'conc': bl_str})

        ## Ref conc
        ref_conc0 = utils.get_value(param.lakes_ref_conc_path, lake_id)
        ref_conc = ref_conc0[indicator]
        mean1 = round(np.exp(ref_conc['log_mean']), 3)
        p5 = round(np.exp(norm.ppf(0.05, loc=ref_conc['log_mean'], scale=ref_conc['log_stdev'])), 3)
        p95 = round(np.exp(norm.ppf(0.95, loc=ref_conc['log_mean'], scale=ref_conc['log_stdev'])), 3)
        stats.append({'name': 'Reference', 'conc': results_str2.format(mean1, p5, p95)})

        ## Current
        stats0 = utils.get_value(param.lakes_wq_stats_path, (lake_id, indicator))

        stats.append({'name': 'Current', 'conc': results_str1.format(round(stats0['median'], 3))})

        ### Scenario calcs
        if trig == 'tbl_calc_btn':
            stats.append({'name': 'Scenario', 'conc': results_str1.format(round(stats0['median']*conc_factor, 3))})

            disabled = False
        else:
            stats.append({'name': 'Scenario', 'conc': 'Press the Run catchment scenario button'})

    else:
        if stats and indicator is not None:
            stats[-1] = {'name': 'Scenario', 'conc': 'Press the Run catchment scenario button'}

        else:
            stats = []

    # print(stats)
    # print(nzsegment)

    return stats, disabled


@callback(
        Output('stats_tbl', 'children'),
        Input('results_tabs', 'value'),
        Input('stats', 'data'),
        State('indicator', 'value'),
        prevent_initial_call=True
        )
def update_stats_table(tab, stats, indicator):

    stats_tbl_data = []
    # print(tab)
    if stats:
        for s in stats:
            if s['name'] in ('Reference', 'Current', 'Scenario'):
                stats_tbl_data.append(s)

    stats_tbl = utils.make_results_table(stats_tbl_data, tab, indicator)

    return stats_tbl


@callback(
        Output('box_plot_fig', 'data'),
        Input('stats', 'data'),
        State('lake_id', 'data'),
        State('indicator', 'value'),
        prevent_initial_call=True
        )
def update_box_plot_fig(stats, lake_id, indicator):

    # print(stats_tbl)
    if stats:
        fig = utils.make_fig_ind(stats, lake_id, indicator)
    else:
        fig = None

    box_plot_fig_enc = utils.encode_obj(fig)

    return box_plot_fig_enc


@callback(
        Output('ts_plot_div', 'children'),
        Input('results_tabs', 'value'),
        Input('box_plot_fig', 'data'),
        prevent_initial_call=True
        )
def update_box_plot(tab, box_plot_fig_enc):

    # print(stats_tbl)
    fig = utils.decode_obj(box_plot_fig_enc)

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
        State('conc_factor', 'data'),
        State('indicator', 'value'),
        State('lake_data', 'data'),
        prevent_initial_call=True,
        )
def make_pdf_report(n_clicks, lc_tbl, stats, lake_id, lake_name, conc_factor, indicator, lake_data):
    """

    """
    ## Pre-calcs
    if indicator == 'ECOLI':
        units = 'CFU/100ml'
    else:
        units = 'mg/l'

    current_scenario = {l['name']: float(l['conc']) for l in stats if l['name'] in ('Current', 'Scenario')}
    improve_perc = int(round((1 - (current_scenario['Scenario']/current_scenario['Current'])) * 100))
    improve_text = 'improvement'
    if improve_perc < 0:
        improve_perc = improve_perc * -1
        improve_text = 'degredation'

    stats0 = utils.get_value(param.lakes_wq_stats_path, (lake_id, indicator))

    obs_count = stats0['count']
    obs_min_date = stats0['from_date']
    obs_max_date = stats0['to_date']
    obs_median = current_scenario['Current']
    scenario_conc = current_scenario['Scenario']
    ref_mean, ref_lower, ref_upper = utils.sep_reference_values([s['conc'] for s in stats if s['name'] == 'Reference'][0])
    nps_check = [True for s in stats if s['name'] == 'Bottom line']

    lc_data = utils.get_value(param.lakes_lc_red_yields_path, lake_id).sort_index()
    lc_data1 = lc_data[[col for col in lc_data.columns if 'reduction' not in col]].round(2).reset_index().copy()
    lc_data1['area_ha'] = lc_data1['area_ha'].round().astype(int)
    tot_area = lc_data1['area_ha'].sum()

    if indicator != 'ECOLI':
        ind_name = param.rivers_indicator_dict[indicator].lower()
    else:
        ind_name = param.rivers_indicator_dict[indicator]
    # river_name = site_data['properties']['Catchment']
    agency = lake_data['properties']['Agency']

    # box_plot_fig = utils.decode_obj(box_plot_fig_enc)
    # print(geometry)
    # geom = orjson.loads(geometry.encode())

    ## plot the catchment polygon and river network
    fig, ax = plt.subplots()

    ## Rivers
    lakes_reaches_gbuf = utils.get_value(param.lakes_catch_reaches_gbuf_path, lake_id)

    lakes_reaches_dict = geobuf.decode(lakes_reaches_gbuf)
    if lakes_reaches_dict is not None:
        rivers_geo = gpd.GeoDataFrame.from_features(lakes_reaches_dict['features'], crs=4326)
        rivers_geo.plot(ax=ax, alpha=0.3)

    ## Catchment
    catch_gbuf = utils.get_value(param.lakes_catch_gbuf_path, lake_id)

    catch_dict = geobuf.decode(catch_gbuf)
    catch_geo = gpd.GeoDataFrame.from_features(catch_dict['features'], crs=4326)
    catch_geo.plot(ax=ax, color='lightgrey', edgecolor='black', alpha=0.6)

    ## Site location
    point = gpd.GeoSeries([Point(lake_data['geometry']['coordinates'])], crs=4326)
    point.plot(ax=ax, color='black')
    # ax.set_xlim(x_min - b, x_max + b)
    # ax.set_ylim(y_min - b, y_max + b)
    # ax.set_extent((x_min - b, y_min - b, x_max + b, y_max + b))
    cx.add_basemap(ax, crs=catch_geo.crs, source=cx.providers.OpenStreetMap['Mapnik'])

    fig.tight_layout()

    ## Make the document
    lc_tbl_dict = {k: l[1] for k, l in param.tbl_cols_dict.items()}

    with tempfile.TemporaryDirectory() as path:
        doc = Document(os.path.join(path, lake_id), geometry_options={"right": "2cm", "left": "2cm"})
        doc.packages.append(Package('array'))
        doc.packages.append(Package('datetime2'))
        doc.packages.append(Package('float'))
        # right_fixed_col = ColumnType('R', 'p{#1}', '\raggedright')
        doc.preamble.append(NoEscape(r'\newcolumntype{R}[1]{>{\raggedleft\arraybackslash}p{#1}}'))
        doc.preamble.append(NoEscape(r'\newcolumntype{L}[1]{>{\raggedright\arraybackslash}p{#1}}'))
        doc.preamble.append(NoEscape(r'\newcolumntype{C}[1]{>{\centering\arraybackslash}p{#1}}'))

        doc.preamble.append(Command('title', f'Scenario builder results for {ind_name} at {lake_name}'))
        # doc.preamble.append(catch_name)
        # doc.preamble.append(Command('date', NoEscape(r'\today')))
        # doc.preamble.append(NoEscape(r'\today'))
        doc.append(NoEscape(r'\maketitle'))

        ## Catchment location section
        with doc.create(Section('Monitoring site and catchment information')):
            text = f'The user-selected site is named "{lake_name}" and monitored by {agency}. The site is also named "{lake_id}" by DOC. '

            # if catch_name != 'Unnamed':
            #     text += f'It is located within the larger {catch_name} catchment '

            doc.append(text)
            doc.append(NoEscape(r'(Figure \ref{fig:catchment}).'))
            # position='htbp'
            with doc.create(Figure(position='H')) as plot:
                plot.add_plot(width=NoEscape(r'0.8\textwidth'))
                plot.add_caption(f'Catchment location with a total area of {tot_area:,} hectares. The black dot is the monitoring site.')
                plot.append(Command('label', 'fig:catchment'))

        ## Section 2
        # sec2 =  Section('Land areas and mitigations')
        sec2 =  Section('Current estimates')
        # sec2.append(NoEscape(r'The current state and user-defined scenario for the land areas and mitigations are shown in Table \ref{tab:mitigations}. The user should refer to the user guide for more details on the background and methodology of the calculations.'))
        doc.append(sec2)

        # 2.1 - Indicator
        sec2a = Subsection('Measured indicator')
        sec_text = f'The user-selected water quality indicator is {ind_name}. The current estimate of the median concentration for the 5 year period from an initial observation at {obs_min_date} to a final observation at {obs_max_date} is {obs_median} {units} '
        sec_text += r'(Table \ref{tab:conc}). '
        sec_text += f'The estimated reference median concentration is {ref_mean} {units} with 5th and 95th percentiles of {ref_lower} and {ref_upper} {units}. The total number of observations during this period was {obs_count}. '
        if nps_check:
            sec_text += r'The relevant NPS-FM 2020 bands associated with the median concentration attribute states are listed in Table \ref{tab:conc}. '
        sec_text += r'Table \ref{tab:mitigations} '
        sec_text += f'lists the land use composition and estimated contribution of each land use to the total load of {ind_name} at the monitoring site. '
        if indicator != 'ECOLI':
            sec_text += r'Table \ref{tab:yields} shows the areas and loss rates of each land use in the selected catchment used in the Scenario Builder WebApp.  '

        sec2a.append(NoEscape(sec_text))

        doc.append(sec2a)

        ## Conc table
        header = ['Name', f'Median concentration {units}']
        base_tbl2 = Table(position="htb")
        base_tbl2.add_caption('Site median concentrations')
        base_tbl2.append(NoEscape(r'\centering'))
        data_table2 = Tabular("| l | R{0.20\linewidth} |")
        data_table2.add_hline()
        data_table2.add_row(header, strict=False)
        data_table2.add_hline()
        for values in stats:
            row = [values['name'], values['conc']]
            data_table2.add_row(row, strict=False)
        data_table2.add_hline()
        base_tbl2.append(data_table2)
        # base_tbl2.add_caption('Site median concentrations')
        base_tbl2.append(Command('label', 'tab:conc'))

        doc.append(base_tbl2)

        ## Land use table
        base_tbl = Table(position="htb")
        base_tbl.add_caption('Land area contributions and mitigations')
        base_tbl.append(NoEscape(r'\centering'))
        data_table = Tabular('| l | R{0.06\linewidth} | R{0.13\linewidth} | R{0.13\linewidth} | R{0.06\linewidth} | R{0.13\linewidth} |')
        data_table.add_hline()
        data_table.add_row(['', MultiColumn(2, align='c|', data='Current'), MultiColumn(3, align='c|', data='Scenario')], strict=False)
        data_table.add_row(list(lc_tbl_dict.values()), strict=False)
        data_table.add_hline()
        # data_table.end_table_header()
        # data_table.add_hline()
        # data_table.add_row((MultiColumn(3, align='r',
        #                     data='Continued on Next Page'),))
        # data_table.add_hline()
        # data_table.end_table_footer()
        # data_table.add_hline()
        # data_table.add_row((MultiColumn(3, align='r',
        #                     data='Not Continued on Next Page'),))
        # data_table.add_hline()
        # data_table.end_table_last_footer()
        for values in lc_tbl:
            row = [str(int(val)) if isinstance(val, (float, int)) else val for val in values.values()]
            data_table.add_row(row, strict=False)
        data_table.add_hline()
        base_tbl.append(data_table)
        # base_tbl.add_caption('Land areas and mitigations')
        base_tbl.append(Command('label', 'tab:mitigations'))

        doc.append(base_tbl)

        ## Land use yields and areas
        if indicator != 'ECOLI':
            base_tbl = Table(position="htb")
            base_tbl.add_caption('Current catchment land use areas and loss rates')
            base_tbl.append(NoEscape(r'\centering'))
            data_table = Tabular('| l | R{0.10\linewidth} | R{0.10\linewidth} | R{0.10\linewidth} |')
            data_table.add_hline()
            data_table.add_row(['', '', MultiColumn(2, align='c|', data='Loss rate (kg/ha/yr)')], strict=False)
            data_table.add_row(['Land use', 'Area (ha)', 'Phosphorus', 'Nitrogen'], strict=False)
            data_table.add_hline()
            # data_table.end_table_header()
            # data_table.add_hline()
            # data_table.add_row((MultiColumn(3, align='r',
            #                     data='Continued on Next Page'),))
            # data_table.add_hline()
            # data_table.end_table_footer()
            # data_table.add_hline()
            # data_table.add_row((MultiColumn(3, align='r',
            #                     data='Not Continued on Next Page'),))
            # data_table.add_hline()
            # data_table.end_table_last_footer()
            for grp, values in lc_data1.iterrows():
                row = values.values.tolist()
                data_table.add_row(row, strict=False)
            data_table.add_hline()
            base_tbl.append(data_table)
            # base_tbl.add_caption('Land areas and mitigations')
            base_tbl.append(Command('label', 'tab:yields'))

            doc.append(base_tbl)

        ## Section 3
        sec3 = Section('Scenario outcomes')
        doc.append(sec3)

        # Indicator
        sec3a = Subsection('Measured indicator')
        sec_text = r'Given the user-defined mitigations/land use change shown in Table \ref{tab:mitigations}, '
        sec_text += f'the estimated {improve_text} at the site is {improve_perc}\% with a resulting median concentration of {scenario_conc} {units} '
        sec_text += r'(Table \ref{tab:conc}). Refer to the user guide for more details on the background and methodology of the calculations.'

        sec3a.append(NoEscape(sec_text))
        doc.append(sec3a)

        doc.generate_pdf(clean_tex=True)

        pdf_path = os.path.join(path, f'{lake_id}.pdf')
        with open(pdf_path, 'rb') as f:
            # pdf_bytes = io.BytesIO(f.read())
            pdf_bytes = f.read()

    return dcc.send_bytes(pdf_bytes, f'{lake_id}.pdf', 'pdf')








