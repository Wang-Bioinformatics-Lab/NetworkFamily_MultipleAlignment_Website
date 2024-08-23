import os
import collections
from typing import List, Tuple, Dict, Set
import json
import uuid
from config import SETS_TEMP_PATH

import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc

from app import app

from alignment import get_data, get_topo_path, get_sets, new_matches, add_pairs, data_for_json


SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["scan", "precursor_mz", "precursor_charge", "mz", "intensity"]
)

# peak tuple for the each peak in a set to store the scan number and peak idx
PeakTuple = collections.namedtuple("PeakTuple", ["scan_num", "peak_idx"])

dash_app = dash.Dash(
    name="setscreation",
    server=app,
    url_base_pathname="/setscreation/",
    external_stylesheets=[dbc.themes.BOOTSTRAP],
)

# dash_app.layout = html.Div([
#     dcc.Location(id='url', refresh=False),
#     dcc.Input(id='task-id-input', type='text', placeholder='Enter task ID', value="c198b31cb3e241ccbf1d7fc2dd9af0c7"),
#     dcc.Input(id='component-input', type='text', placeholder='Enter component number', value="1"),
#     html.Button('Process and Save JSON', id='process-button'),
#     html.Div(id='output-path')
# ])


dash_app.layout = dbc.Container(
    [
        dbc.Row(
            dbc.Col(
                html.H1("Molecular Networking Peak Alignment", style={'textAlign': 'center',  'padding': '10px'}),
                width=12
            )
        ),
        dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Row(html.H3("Data Selection", className="card-title", style={'fontSize': '18px'}))
                ),
                dbc.CardBody(
                    [
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        dbc.Label("Task ID", html_for="task-id-input"),
                                        dbc.Input(
                                            id='task-id-input',
                                            type='text',
                                            placeholder='Enter task ID',
                                            value="c198b31cb3e241ccbf1d7fc2dd9af0c7",
                                            # style={'width': '100%'},
                                            size='sm'
                                        ),
                                    ],
                                    width=6
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label("Component Number", html_for="component-input"),
                                        dbc.Input(
                                            id='component-input',
                                            type='text',
                                            placeholder='Enter component number',
                                            value="1",
                                            # className="mb-4",
                                            # style={'fontSize': '16px'},
                                            size="sm" 
                                        ),
                                    ],
                                    width=4,
                                ),
                                dbc.Col(
                                    dbc.Button(
                                        'Process and Save JSON',
                                        id='process-button',
                                        color='primary',
                                        size="sm"
                                    ),
                                    width=2,
                                    className="d-flex align-items-end"
                                ),
                            ],
                            style={'margin-bottom': '10px'}
                        ),
                        dbc.Row(
                            dbc.Col(
                                html.Div(id='output-path'),
                                width=12
                            )
                        )
                    ]
                )
            ]
        ),
    ],
    fluid=True
)

@dash_app.callback(
    Output('output-path', 'children'),
    [Input('process-button', 'n_clicks')],
    [State('task-id-input', 'value'), State('component-input', 'value')]
)
def process_and_save_json(n_clicks, task_id, component_number):
    if n_clicks is None or not task_id or not component_number:
        return dash.no_update

    # Figuring out the full path hashed on task_id and component number
    SAVE_PATH = os.path.join(SETS_TEMP_PATH, f"{task_id}_{component_number}.json")
    
    filtered_spec_dic, df_comp, scans = get_data(int(component_number), task_id)
    topo_path, alignments = get_topo_path(filtered_spec_dic, df_comp, scans)
    transitive_sets = get_sets(topo_path, alignments)
    new_pairs = new_matches(topo_path, filtered_spec_dic)
    final_sets = add_pairs(transitive_sets, new_pairs)
    json_sets, new_spec_dic = data_for_json(final_sets, filtered_spec_dic, topo_path)

    with open(SAVE_PATH, 'w') as f:
        # json.dump((component, json_sets, new_spec_dic), f)
        json.dump((int(component_number), json_sets, new_spec_dic), f)
    
    output_result = []


    # Create a linkout to the other page with the json file in the url
    # linkout = dash.dcc.Link('View Alignment', href=f'/spectraalignment?filename={task_id}_{component_number}.json', target='_blank')

    # create a link button
    # linkout = dbc.Button('View Alignment', href=f'/spectraalignment?filename={task_id}_{component_number}.json', target='_blank', color='primary', size = 'sm')
    linkout = html.A(
        dbc.Button('View Alignment', color='primary', size = 'sm'),
        href=f'/spectraalignment?filename={task_id}_{component_number}.json',
        target='_blank'
    )

    output_result.append(html.P(f"Saved to {SAVE_PATH}"))
    output_result.append(linkout)

    return output_result
    

# Setting file-name-input value from the url search parameters
@dash_app.callback(
    [
        Output('task-id-input', 'value'),
        Output('component-input', 'value')
     ],
    Input('url', 'search')
)
def update_file_name(search):

    try:
        # Parsing out the search field to grab the file name
        import urllib.parse
        params_dict = urllib.parse.parse_qs(search[1:])
        return [params_dict['task'][0], params_dict['component'][0]]
    except:
        return [dash.no_update, dash.no_update]
    
if __name__ == '__main__':
    app.run(debug=True)