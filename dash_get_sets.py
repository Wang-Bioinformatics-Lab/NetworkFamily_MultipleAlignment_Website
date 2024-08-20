import requests
import pandas as pd
import io
import os
import collections
from typing import List, Tuple, Dict, Set
import json
import uuid

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

dash_app.layout = html.Div([
    dcc.Input(id='task-id-input', type='text', placeholder='Enter task ID', value="c198b31cb3e241ccbf1d7fc2dd9af0c7"),
    dcc.Input(id='component-input', type='text', placeholder='Enter component number', value="1"),
    html.Button('Process and Save JSON', id='process-button'),
    html.Div(id='output-path')
])

SETS_TEMP_PATH = "./temp/sets"

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
    linkout = dash.dcc.Link('View Alignment', href=f'/spectraalignment?filename={task_id}_{component_number}.json', target='_blank')

    output_result.append(html.P(f"Saved to {SAVE_PATH}"))
    output_result.append(linkout)

    return output_result
    
if __name__ == '__main__':
    app.run(debug=True)