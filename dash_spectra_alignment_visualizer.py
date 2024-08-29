# -*- coding: utf-8 -*-
import dash
from dash import dcc
from dash import html
from dash import Dash, callback, Output, Input, State, dash_table

import dash_bootstrap_components as dbc

import plotly.graph_objects as go
from dash.dependencies import Input, Output, State, ALL
import collections
import json
import os

from config import SETS_TEMP_PATH
import numpy as np

SpectrumTuple = collections.namedtuple('SpectrumTuple', ['scan', 'precursor_mz', 'precursor_charge', 'mz', 'intensity'])
PeakTuple = collections.namedtuple('PeakTuple', ['scan_num', 'peak_idx'])

from app import app
dash_app = dash.Dash(
    name="spectraalignment",
    server=app,
    url_base_pathname="/spectraalignment/",
    external_stylesheets=[dbc.themes.BOOTSTRAP],
)

dash_app.title = "Spectral Alignment Visualization"

# setting tracking token
dash_app.index_string = """<!DOCTYPE html>
<html>
    <head>
        <!-- Umami Analytics -->
        <script async defer data-website-id="ENTER YOUR TOKEN HERE" src="https://analytics.gnps2.org/umami.js"></script>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>"""

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


def calculate_percent_top_10(peak_sets, spec_dic):
    top_10_dic = {}

    # change spectrum tuple to dictionary
    for key, s in spec_dic.items():
        # get top 10 peaks based on intensity
        top_10_indices = np.argsort(s.intensity)[-10:][::-1]
        top_10_mz = [s.mz[i] for i in top_10_indices]
        top_10_intensity = [s.intensity[i] for i in top_10_indices]

        # use a set of tuples (mz, intensity)
        top_10_peaks = set(zip(top_10_mz, top_10_intensity))

        top_10_dic[int(s.scan)] = {
            "scan": int(s.scan),
            "top_10_peaks": top_10_peaks
        }

    percent_top_10_info = []

    # iterate over each set
    for set_index, match_list in enumerate(peak_sets):
        top_10_count = 0  # count for peaks that are in the top 10
        total_peaks = len(match_list) 
        
        # iterate over each peak in the current set
        for peak_tuple in match_list:
            spectrum_id = peak_tuple.scan_num  # get scan num
            peak_index = peak_tuple.peak_idx  # get peak index 
            
            # get top 10 peaks for this spectrum
            top_10_peaks = top_10_dic[spectrum_id]['top_10_peaks']
            
            # get mz and intensity of the current peak
            peak_mz = spec_dic[int(spectrum_id)].mz[int(peak_index)]
            peak_intensity = spec_dic[spectrum_id].intensity[peak_index]

            # check if the peak is in the top 10
            if (peak_mz, peak_intensity) in top_10_peaks:
                top_10_count += 1

        # calculate the percentage of peaks in the set that are in the top 10
        if total_peaks > 0:
            percent_top_10 = (top_10_count / total_peaks) * 100
        else:
            percent_top_10 = 0
        
    # add the result to the list
        percent_top_10_info.append({
            'set_number': set_index + 1,
            'set_size': total_peaks,
            'percent': percent_top_10
        })

    # create a table
    table = dash_table.DataTable(
        data=percent_top_10_info,
        columns=[
            {'name': 'Set Number', 'id': 'set_number', 'type': 'numeric'},
            {'name': 'Set Size', 'id': 'set_size', 'type': 'numeric'},
            {'name': 'Percent', 'id': 'percent', 'type': 'numeric', 'format': {'specifier': '.2f'}},
        ],
        sort_action='native',  # allow for sorting
        style_table={'overflowX': 'auto'},
        style_header={'fontWeight': 'bold'},
        style_cell={'textAlign': 'left', 'padding': '5px'},
        style_as_list_view=True,
        page_size=10  # 10 rows per page
    )

    return table


dash_app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.H1('Molecular Networking Peak Alignment', style={'textAlign': 'center',  'padding': '10px'}),

    # dbc card for data input/custom order/selection dropdown
    dbc.Card(
        [
            # header section
            dbc.CardHeader(
                dbc.Row([
                    dbc.Col(html.H3("Data Input", className="card-title", style={'fontSize': '20px'})),
                    dbc.Col(dbc.Button("Show/Hide", id="toggle-button", color="secondary", size="sm"), width="auto")
                ], align="center"),
                style={'border-bottom': '1px solid #ccc'} 
            ),

            # body
            dbc.Collapse(
                dbc.CardBody(
                    dbc.Row([
                        dbc.Col([
                            dbc.Label("File Name", html_for='file-name-input'),
                            dcc.Input(id='file-name-input', type='text', placeholder='Enter the file name', 
                                      style={'width': '100%', 'fontSize': '14px', 'padding': '5px', 'borderRadius': '5px', 'border': '1px solid #ccc'}),
                        ], width=12, style={'margin-bottom': '10px'}),

                        dbc.Col([
                            dbc.Label("Custom Order of Scan Numbers (Optional)", html_for='custom-order-input'),
                            dcc.Input(id='custom-order-input', type='text', placeholder='Enter custom order of scan numbers separated by commas', 
                                      style={'width': '100%', 'fontSize': '14px', 'padding': '5px', 'borderRadius': '5px', 'border': '1px solid #ccc'}),
                        ], width=12, style={'margin-bottom': '10px'}),

                        dbc.Col([
                            dbc.Label("Select Sorting Order (Optional)", html_for='sort-order-dropdown'),
                            dcc.Dropdown(
                                id='sort-order-dropdown',
                                options=[
                                    {'label': 'Default (Topological Sort) Order', 'value': 'original'},
                                    {'label': 'Ascending by Precursor m/z', 'value': 'asc'},
                                    {'label': 'Descending by Precursor m/z', 'value': 'desc'},
                                    {'label': 'Custom Order', 'value': 'custom'}
                                ],
                                placeholder='Select sorting order', 
                                style={'fontSize': '14px'}
                            ),
                        ], width=12, style={'margin-bottom': '10px'}),

                        dbc.Col([
                            dbc.Label("X-axis Lower Limit (Optional)", html_for='lower-lim-input'),
                            dcc.Input(id='lower-lim-input', type='text', placeholder='Enter the minimum mz', 
                                      style={'width': '100%', 'fontSize': '14px', 'padding': '5px', 'borderRadius': '5px', 'border': '1px solid #ccc'})
                        ], width=12, style={'margin-bottom': '10px'}),

                        dbc.Col([
                            dbc.Label("X-axis Upper Limit (Optional)", html_for='upper-lim-input'),
                            dcc.Input(id='upper-lim-input', type='text', placeholder='Enter the maximum mz', 
                                      style={'width': '100%', 'fontSize': '14px', 'padding': '5px', 'borderRadius': '5px', 'border': '1px solid #ccc'})
                        ], width=12, style={'margin-bottom': '10px'})
                    ])
                ),
                id="collapse-input",
                is_open=True,
            ),
        ],
        style={'margin-bottom': '20px'}
    ),

    # card for set top 10 percents
    dbc.Card(
        [
            dbc.CardHeader(
                dbc.Row([
                    dbc.Col(html.H3("Set Info with Top 10 Peak Intensity Percentages", className="card-title", style={'fontSize': '20px'})),
                    dbc.Col(dbc.Button("Show/Hide", id="toggle-set-percents", color="secondary", size="sm"), width="auto")
                ], align="center"),
                style={'border-bottom': '1px solid #ccc'}  
            ),

            dbc.Collapse(
                dbc.CardBody(
                    dbc.Row([
                        html.Div(id='set_info', style={'margin-bottom': '10px'}), 
                        html.Div(id='percent-top-10-output')
                    ])
                    
                ),
                id="collapse-set-percents",
                is_open=False,
            ),
        ],
        style={'margin-bottom': '20px'}
    ),
    
    # card for set info
    dbc.Card(
        [
            dbc.CardHeader(
                dbc.Row([
                    dbc.Col(html.H3("Set Info", className="card-title", style={'fontSize': '20px'})),
                    dbc.Col(dbc.Button("Show/Hide", id="toggle-set-info", color="secondary", size="sm"), width="auto")
                ], align="center"),
                style={'border-bottom': '1px solid #ccc'}  
            ),

            dbc.Collapse(
                dbc.CardBody(
                    html.Div(id='set-info')
                ),
                id="collapse-set-info",
                is_open=False,
            ),
        ],
        style={'margin-bottom': '20px'}
    ),

    dbc.Col(dbc.Button('Display Spectra', id='display-button', n_clicks=0, color="primary", className="ml-3"), width="auto", style = {'margin-bottom': '20px'}),
    html.Div(id='graphs-container', style={'padding': '0', 'margin': '0'}),
    dcc.Store(id='clicked-peak-store', data={'scan': None, 'peak_idx': None}),
    dcc.Store(id='highlighted-sets', data=[]),
    dcc.Store(id='spec-data', data=None),
    html.Div(style={'height': '100px'})
])

# function to create a spectrum figure
def make_spectrum_fig(spectrum, spec_id, highlighted_sets, clicked_peak, show_x_labels, max_mz, sets, spec_dic, lower_lim, upper_lim):
    fig = go.Figure()
    mz = spectrum.mz
    intensities = spectrum.intensity

    # initialize all peak colors to grey
    colors = ['grey'] * len(mz)
    annotations = []

    # highlight top 3 largest sets in red, green, and blue
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

        # highlight matched peaks, red is common (marked as 2), blue is shifted (marked as 1)
        for match_set in sets:
            if (clicked_scan, clicked_idx) in match_set:
                temp_dic = {}
                for peak in match_set:
                    # if it's an exact match
                    if spec_dic[peak[0]].mz[peak[1]] in temp_dic:
                        temp_dic[spec_dic[peak[0]].mz[peak[1]]] = 2
                    # if the mz is within the tolerance of 0.1
                    elif (spec_dic[peak[0]].mz[peak[1]] + 0.1) in temp_dic:
                        temp_dic[spec_dic[peak[0]].mz[peak[1]]] = 2
                    elif (spec_dic[peak[0]].mz[peak[1]] - 0.1) in temp_dic:
                        temp_dic[spec_dic[peak[0]].mz[peak[1]]] = 2
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
        margin=dict(t=0, b=0, l=0, r=20)
    )

    # set range for x-axis 
    if lower_lim and upper_lim:
        fig.update_xaxes(range=[lower_lim, upper_lim], showticklabels=show_x_labels)
    elif lower_lim:
        fig.update_xaxes(range=[lower_lim, max_mz], showticklabels=show_x_labels)
    elif upper_lim:
        fig.update_xaxes(range=[0, upper_lim], showticklabels=show_x_labels)
    else:
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
    # Cleaning up the filename
    file_name = os.path.join(SETS_TEMP_PATH, file_name)

    # loading spectra
    peak_sets, spec_dic, max_mz, max_size = _load_peaksets(file_name)

    custom_order_list = [int(scan) for scan in custom_order.split(',')]
    ordered_spec_dic = {scan: spec_dic[scan] for scan in custom_order_list if scan in spec_dic}

    for i, data in enumerate(clickData):
        if data and 'points' in data:
            point = data['points'][0]
            spec_id = list(ordered_spec_dic.keys())[i]
            peak_idx = point['pointIndex']
            return {'scan': spec_id, 'peak_idx': peak_idx}
    return current_data 

@dash_app.callback(
    Output('set-info', 'children'),
    Input('clicked-peak-store', 'data'),
    State('file-name-input', 'value'),
    State('custom-order-input', 'value'),
    prevent_initial_call=True
)
def display_set_info(clicked_peak, file_name, custom_order):
    if clicked_peak['scan'] is None or clicked_peak['peak_idx'] is None:
        return "Click on a peak to see its set information."
    
    file_name = os.path.join(SETS_TEMP_PATH, file_name)

    # load data again 
    peak_sets, spec_dic, max_mz, max_size = _load_peaksets(file_name)

    clicked_scan = clicked_peak['scan']
    clicked_idx = clicked_peak['peak_idx']

    # find the set containing the clicked peak
    peak_set = None
    for i, p_set in enumerate(peak_sets):
        if (clicked_scan, clicked_idx) in p_set:
            peak_set = p_set
            set_num = i + 1
            break

    if peak_set is None:
        return "No matching set found for the clicked peak."
    
    # get custom order of scan numbers
    custom_order_list = [int(scan) for scan in custom_order.split(',')]

    # make dictionary of peaks by scan number
    peaks_dict = {peak[0]: peak for peak in peak_set}

    # get peaks in the custom order
    ordered_peaks = [peaks_dict[scan] for scan in custom_order_list if scan in peaks_dict]

    data = [
        {
            'Scan #': peak[0],
            'Peak Index': peak[1],
            'Precursor m/z': f"{spec_dic[peak[0]].precursor_mz:.3f}",
            'Peak m/z': f"{spec_dic[peak[0]].mz[peak[1]]:.3f}",
            'Peak Intensity': f"{spec_dic[peak[0]].intensity[peak[1]]:.3f}", 
            'Set Number': set_num
        }
        for peak in ordered_peaks
    ]

    # create table
    table = dash_table.DataTable(
        data = data,
        columns = [
            {'name': 'Scan #', 'id': 'Scan #', 'type': 'numeric'},
            {'name': 'Peak Index', 'id': 'Peak Index', 'type': 'numeric'},
            {'name': 'Precursor m/z', 'id': 'Precursor m/z', 'type': 'numeric'},
            {'name': 'Peak m/z', 'id': 'Peak m/z', 'type': 'numeric'},
            {'name': 'Peak Intensity', 'id': 'Peak Intensity', 'type': 'numeric'},
            {'name': 'Set Number', 'id': 'Set Number', 'type': 'numeric'},
        ],
        sort_action='native',
        style_table={'overflowX': 'auto'},
        style_header={'fontWeight': 'bold'},
        style_cell={'textAlign': 'left', 'padding': '5px'},
        style_as_list_view=True,
        page_size=10 # 10 rows per page
    )

    return table

# callback for displaying spectra
@dash_app.callback(
    [Output('percent-top-10-output', 'children'),
     Output('set_info', 'children'), 
     Output('graphs-container', 'children'), 
     Output('highlighted-sets', 'data'),
     Output('custom-order-input', 'value')],
    Input('display-button', 'n_clicks'),
    Input('clicked-peak-store', 'data'),
    State('file-name-input', 'value'),
    State('custom-order-input', 'value'),
    State('sort-order-dropdown', 'value'),
    Input('lower-lim-input', 'value'),
    Input('upper-lim-input', 'value'),
    prevent_initial_call=True
)

def display_spectra(n_clicks, clicked_peak, file_name, custom_order, sort_order, lower_lim, upper_lim):    
    file_name = os.path.join(SETS_TEMP_PATH, file_name)

    peak_sets, spec_dic, max_mz, max_size = _load_peaksets(file_name)

    # get largest sets
    largest_sets = sorted(peak_sets, key=len, reverse=True)[:3]

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
        fig = make_spectrum_fig(spectrum, scan, largest_sets, clicked_peak, show_x_labels, max_mz, peak_sets, spec_dic, lower_lim, upper_lim)
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
                        'text-align': 'right', 
                        'width': '230px' # adjust if needed 
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

    # calculate the percent top 10 for the sets
    percent_top_10_info = calculate_percent_top_10(peak_sets, spec_dic)

    sets_info = f'The table below displays the information for {len(peak_sets)} sets, as well as the percentage of peaks within each set that are in the top 10 intensities per spectrum'

    return percent_top_10_info, sets_info, graphs, largest_sets, current_order_str

# Setting file-name-input value from the url search parameters
@dash_app.callback(
    Output('file-name-input', 'value'),
    Input('url', 'search')
)
def update_file_name(search):

    try:
        # Parsing out the search field to grab the file name
        import urllib.parse
        params_dict = urllib.parse.parse_qs(search[1:])
        return params_dict['filename'][0]
    except:
        return ""

# callback to toggle collapse for data input
@dash_app.callback(
    Output("collapse-input", "is_open"),
    Input("toggle-button", "n_clicks"),
    State("collapse-input", "is_open"),
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open

# call back for collapsing largest sets 
@dash_app.callback(
    Output("collapse-largest-sets", "is_open"),
    Input("toggle-largest-sets", "n_clicks"),
    State("collapse-largest-sets", "is_open")
)
def toggle_largest_sets(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open

# call back for collapsing set info
@dash_app.callback(
    Output("collapse-set-info", "is_open"),
    Input("toggle-set-info", "n_clicks"),
    State("collapse-set-info", "is_open")
)
def toggle_largest_sets(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open

# call back for collapsing set percents
@dash_app.callback(
    Output("collapse-set-percents", "is_open"),
    Input("toggle-set-percents", "n_clicks"),
    State("collapse-set-percents", "is_open")
)
def toggle_largest_sets(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open

if __name__ == '__main__':
    app.run(debug=True)