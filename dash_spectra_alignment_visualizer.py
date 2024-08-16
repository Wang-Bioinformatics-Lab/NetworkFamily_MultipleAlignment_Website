# -*- coding: utf-8 -*-
import dash
from dash import dcc
from dash import html
from dash import Dash, html, dcc, callback, Output, Input, State

import dash_bootstrap_components as dbc

import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output, State, ALL
from plotly.subplots import make_subplots
from typing import List, Tuple, Dict, Set
import collections
import json

SpectrumTuple = collections.namedtuple('SpectrumTuple', ['scan', 'precursor_mz', 'precursor_charge', 'mz', 'intensity'])
PeakTuple = collections.namedtuple('PeakTuple', ['scan_num', 'peak_idx'])

from app import app
dash_app = dash.Dash(
    name="spectraalignment",
    server=app,
    url_base_pathname="/spectraalignment/",
    external_stylesheets=[dbc.themes.BOOTSTRAP],
)

dash_app.title = "Dash Interface 1"

# turn tuple back into peak tuple
def tuple_to_peak_tuple(tup):
    return PeakTuple(*tup)

def _load_peaksets(file_name):
    # read in data from json file
    with open(file_name, 'r') as openfile:
        component, json_sets, new_spec_dic = json.load(openfile)

    # add peak tuples to sets
    peak_sets = []
    for json_match_list in json_sets:
        peak_sets.append([tuple_to_peak_tuple(tup) for tup in json_match_list])

    spec_dic = {}
    # change spectra to dictionary format
    for key, s in new_spec_dic.items():
        spec_dic[int(key)] = SpectrumTuple(s['scan'], s['precursor_mz'], s['precursor_charge'], s['mz'], s['intensity'])

    # get largest mz value to set limits for x-axis
    max_mz = max(max(spec.mz) for spec in spec_dic.values())

    # get largest size of set
    max_size = max(len(s) for s in peak_sets)

    return peak_sets, spec_dic, max_mz, max_size

dash_app.layout = html.Div([
    html.H1('Molecular Networking Peak Alignment', style={'textAlign': 'center'}),
    
    dcc.Input(id='file-name-input', type='text', placeholder='Enter the file name', style={'width': '50%'}),
    html.Button('Display Spectra', id='display-button', n_clicks=0),

    dcc.Input(id='custom-order-input', type='text', placeholder='Enter custom order of scan numbers separated by commas', style={'width': '50%'}),
    dcc.Dropdown(
        id='sort-order-dropdown',
        options=[
            {'label': 'Default (Topological Sort) Order', 'value': 'original'},
            {'label': 'Ascending by Precursor m/z', 'value': 'asc'},
            {'label': 'Descending by Precursor m/z', 'value': 'desc'},
            {'label': 'Custom Order', 'value': 'custom'}
        ],
        placeholder='Select sorting order'
    ),
    html.Div(id='file-info'),
    html.Div(id='largest-sets'),
    html.Div(id='graphs-container', style={'padding': '0', 'margin': '0'}),
    dcc.Store(id='clicked-peak-store', data={'scan': None, 'peak_idx': None}),
    dcc.Store(id='highlighted-sets', data=[]),
    dcc.Store(id='spec-data', data=None),
    html.Div(style={'height': '100px'})
])

# function to create a spectrum figure
def make_spectrum_fig(spectrum, spec_id, highlighted_sets, clicked_peak, show_x_labels, max_mz, sets, spec_dic):
    fig = go.Figure()
    mz = spectrum.mz
    intensities = spectrum.intensity

    # initialize all peak colors to grey
    colors = ['grey'] * len(mz)
    annotations = []

    # jighlight top 3 largest sets in red, green, and blue
    for idx, highlight_set in enumerate(highlighted_sets):
        color = ['red', 'green', 'blue'][idx % 3]
        for scan_num, peak_idx in highlight_set:
            if scan_num == spec_id:
                colors[peak_idx] = color

    # if a peak is clicked, change its color and the color of matched peaks to red
    if clicked_peak['scan'] is not None and clicked_peak['peak_idx'] is not None:
        clicked_scan = clicked_peak['scan']
        clicked_idx = clicked_peak['peak_idx']

        # change all peaks to grey again
        colors = ['grey'] * len(mz)

        # highlight the clicked peak
        if clicked_scan == spec_id:
            colors[clicked_idx] = 'red'
            annotations.append(
                dict(
                    x=mz[clicked_idx],
                    y=intensities[clicked_idx],
                    text=f'm/z: {float(mz[clicked_idx]):.3f}, Intensity: {float(intensities[clicked_idx]):.3f}',
                    showarrow=True,
                    arrowhead=2
                )
            )

        # highlight matched peaks
        for match_set in sets:
            if (clicked_scan, clicked_idx) in match_set:
                temp_dic = {}
                for peak in match_set:
                    if spec_dic[peak[0]].mz[peak[1]] in temp_dic:
                        temp_dic[spec_dic[peak[0]].mz[peak[1]]] += 1
                    else:
                        temp_dic[spec_dic[peak[0]].mz[peak[1]]] = 1
                for peak in match_set:
                    if peak[0] == spec_id:
                        colors[peak[1]] = 'blue' if temp_dic[spec_dic[peak[0]].mz[peak[1]]] == 1 else 'red'
                        annotations.append(
                            dict(
                                x=mz[peak[1]],
                                y=intensities[peak[1]],
                                text=f'm/z: {float(mz[peak[1]]):.3f}, Intensity: {float(intensities[peak[1]]):.3f}',
                                showarrow=True,
                                arrowhead=2
                            )
                        )

    # plot the peaks as a bar chart
    fig.add_trace(
        go.Bar(
            x=mz,
            y=intensities,
            name=f'Spectrum {spec_id}',
            hoverinfo='x+y',
            marker_color=colors,
            width=1
        )
    )

    fig.update_layout(
        annotations=annotations,
        yaxis_title='Intensity',
        barmode='overlay',
        showlegend=False,
        height=100,  # adjust if needed
        margin=dict(t=0, b=0, l=0, r=0)
    )

    # set range for x-axis to make all graphs consistent
    fig.update_xaxes(range=[0, max_mz], showticklabels=show_x_labels)

    return fig

# callback for clicked peak
@dash_app.callback(
    Output('clicked-peak-store', 'data'),
    Input({'type': 'spectrum-bar', 'index': ALL}, 'clickData'),
    State('clicked-peak-store', 'data'),
    State('file-name-input', 'value'),
    State('custom-order-input', 'value'),
    prevent_initial_call=True
)
def update_clicked_peak(clickData, current_data, file_name, custom_order):
    # loading spectra
    peak_sets, spec_dic, max_mz, max_size = _load_peaksets(file_name)

    custom_order_list = [int(scan) for scan in custom_order.split(',')]
    ordered_spec_dic = {scan: spec_dic[scan] for scan in custom_order_list if scan in spec_dic}

    for i, data in enumerate(clickData):
        if data and 'points' in data:
            point = data['points'][0]
            # spec_id = list(spec_dic.keys())[i]
            spec_id = list(ordered_spec_dic.keys())[i]
            peak_idx = point['pointIndex']
            return {'scan': spec_id, 'peak_idx': peak_idx}
    return current_data

# callback for displaying spectra
@dash_app.callback(
    [Output('largest-sets', 'children'), 
     Output('graphs-container', 'children'), 
     Output('highlighted-sets', 'data'),
     Output('custom-order-input', 'value')],
    Input('display-button', 'n_clicks'),
    Input('clicked-peak-store', 'data'),
    State('file-name-input', 'value'),
    State('custom-order-input', 'value'),
    State('sort-order-dropdown', 'value'),
    prevent_initial_call=True
)

def display_spectra(n_clicks, clicked_peak, file_name, custom_order, sort_order):    
    peak_sets, spec_dic, max_mz, max_size = _load_peaksets(file_name)

    # get largest sets
    largest_sets = sorted(peak_sets, key=len, reverse=True)[:3]

    # get largest sets information
    largest_sets_info = []
    for i, s in enumerate(largest_sets, 1):
        peaks_info = ', '.join([f'Scan: {p[0]} Peak Index: {p[1]}' for p in s])
        largest_sets_info.append(html.Div(f'Largest Set {i}, size {len(s)}: {peaks_info}'))

    ordered_spec_dic = spec_dic

    if sort_order and sort_order == 'original':
        ordered_spec_dic = spec_dic
    elif sort_order == 'desc':
        ordered_spec_dic = dict(sorted(ordered_spec_dic.items(), key=lambda item: item[1].precursor_mz, reverse=(sort_order == 'desc')))
    elif sort_order == 'asc':    
        ordered_spec_dic = dict(sorted(ordered_spec_dic.items(), key=lambda item: item[1].precursor_mz, reverse=(sort_order == 'desc')))
    elif sort_order == 'custom':
        # get custom order
        if custom_order:
            try:
                custom_order_list = [int(scan) for scan in custom_order.split(',')]
            except ValueError:
                custom_order_list = []
        else:
            custom_order_list = []
        ordered_spec_dic = {scan: spec_dic[scan] for scan in custom_order_list if scan in spec_dic}
    else:
        ordered_spec_dic = spec_dic

    graphs = []

    for i, (scan, spectrum) in enumerate(ordered_spec_dic.items()):
        # only show x axis label for the last spectrum
        show_x_labels = (i == len(spec_dic) - 1)
        fig = make_spectrum_fig(spectrum, scan, largest_sets, clicked_peak, show_x_labels, max_mz, peak_sets, spec_dic)
        graphs.append(
            html.Div([
                html.Div(
                    dcc.Markdown(
                    f'Scan {scan}  Precur m/z: {spec_dic[scan].precursor_mz:.3f}'),
                    style={
                        'transform': 'rotate(0deg)',
                        'height': '100px',
                        'margin-right': '10px',
                        'white-space': 'nowrap',
                        'text-align': 'right',  # adjust if needed
                        'width': '230px'
                    }
                ),
                dcc.Graph(
                    figure=fig,
                    id={'type': 'spectrum-bar', 'index': scan},
                    style={'flex': '1'}
                )
            ], style={'display': 'flex', 'align-items': 'center', 'margin-bottom': '0px'})
        )

    # get scan numbers in the current order
    current_order = list(ordered_spec_dic.keys())
    current_order_str = ', '.join(map(str, current_order))

    return largest_sets_info, graphs, largest_sets, current_order_str

# Setting file-name-input value from the url search parameters
@dash_app.callback(
    Output('file-name-input', 'value'),
    Input('url', 'search')
)
def update_file_name(search):
    return search[1:] if search else ''

    # Parsing out the search field to grab the file name
    return "XXX"



if __name__ == '__main__':
    app.run(debug=True)