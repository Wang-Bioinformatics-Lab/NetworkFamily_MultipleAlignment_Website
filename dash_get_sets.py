import requests
import pandas as pd
import io
import os
from matchms.importing import load_from_mgf
from matchms import calculate_scores
from matchms.similarity import ModifiedCosine
import networkx as nx
from matchms.filtering import normalize_intensities, remove_peaks_around_precursor_mz, remove_peaks_outside_top_k
import collections
import numpy as np
from typing import List, Tuple, Dict, Set
import json
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State

def get_data(component, task_id):
    pairs_url = f'https://gnps2.org/resultfile?task={task_id}&file=nf_output/networking/pairs_with_components.tsv'
    pairs_response = requests.get(pairs_url)
    if pairs_response.status_code != 200:
        return 'Failed to retrieve TSV data'
    
    pairs_data = io.BytesIO(pairs_response.content)
    df = pd.read_csv(pairs_data, sep='\t')
    # print(df.columns)

    # get MGF file
    mgf_url = f'https://gnps2.org/result?task={task_id}&viewname=specms&resultdisplay_type=task'
    mgf_response = requests.get(mgf_url)
    if mgf_response.status_code != 200:
        return 'Failed to retrieve MGF data'

    mgf_data_str = mgf_response.content.decode('utf-8')
    mgf_data = io.StringIO(mgf_data_str)
    spectrums = list(load_from_mgf(mgf_data))
    
    # process DataFrame
    df_comp = df[df['ComponentIndex'] == component]
    # print(df_comp)
    df_comp = df_comp[['CLUSTERID1', 'CLUSTERID2', 'ComponentIndex', 'Cosine']]
    scan_numbers1 = set(df_comp['CLUSTERID1'].astype(int).tolist())
    scan_numbers2 = set(df_comp['CLUSTERID2'].astype(int).tolist())
    scan_numbers = set(scan_numbers2) | set(scan_numbers1)

    # get spectra that match the scan numbers
    filtered_spectra = [s for s in spectrums if int(s.metadata.get("scans")) in scan_numbers]

    # preprocess data
    preprocessed_spectra = []
    for spectrum in filtered_spectra:
        spectrum = normalize_intensities(spectrum)
        spectrum = remove_peaks_around_precursor_mz(spectrum_in = spectrum, mz_tolerance = 17)
        spectrum = remove_peaks_outside_top_k(spectrum_in = spectrum, k = 6, mz_window = 25)
        preprocessed_spectra.append(spectrum)
    
    return preprocessed_spectra

def get_final_sets(component, preprocessed_spectra):
    # calculate modifed cosine scores
    similarity_measure = ModifiedCosine(tolerance = 0.1, mz_power = 0.0, intensity_power = 1.0) 
    scores = calculate_scores(preprocessed_spectra, preprocessed_spectra, similarity_measure)

    # print(f"(Scores names: {scores.score_names})")

    # make an undirected graph for mst
    graph = nx.Graph()
    # add nodes/spectra
    for i, spectrum in enumerate(preprocessed_spectra):
        # add index i and scan number for the node
        graph.add_node(i, scan = int(spectrum.metadata.get("scans")))

    # add edges/modified cosine similarities
    for (spec1, spec2, score) in scores:
        if spec1 != spec2:
            graph.add_edge(spec1, spec2, weight = 1 - score[0])
            
    # kruskal's min spanning tree algorithm
    mst = nx.minimum_spanning_tree(graph)

    # output mst edges
    # print("MST Edges:")
    # for edge in mst.edges(data=True):
    #     spec1, spec2, weight = edge
    #     print(f"(Scan {spec1.get('scans')}) - "
    #         f"(Scan {spec2.get('scans')}), "
    #         f"Similarity Score: {1 - weight['weight']}") # 1 - weight to get the actual cosine score
        
    # create spectrum tuple to store them in a dictionary
    SpectrumTuple = collections.namedtuple("SpectrumTuple", ["scan", "precursor_mz", "precursor_charge", "mz", "intensity"])
    spec_dic = {}

    # change spectra to dictionary format
    for s in preprocessed_spectra:
        spec_dic[int(s.metadata.get('scans'))] = SpectrumTuple(s.metadata.get('scans'), s.metadata.get('precursor_mz'), 
                                                        s.metadata.get('charge'), 
                                                        s.peaks.mz, 
                                                        s.peaks.intensities)
        

    # xianghu's _cosine_fast function to get cosine scores and peak matches
    def _cosine_fast(
            spec: SpectrumTuple,
            spec_other: SpectrumTuple,
            fragment_mz_tolerance: float,
            allow_shift: bool,
    ) -> Tuple[float, List[Tuple[int, int]]]:
        precursor_charge = max(spec.precursor_charge, 1)
        precursor_mass_diff = (
                                    spec.precursor_mz - spec_other.precursor_mz
                            ) * precursor_charge
        # Only take peak shifts into account if the mass difference is relevant.
        num_shifts = 1
        if allow_shift and abs(precursor_mass_diff) >= fragment_mz_tolerance:
            num_shifts += precursor_charge
        other_peak_index = np.zeros(num_shifts, np.uint16)
        mass_diff = np.zeros(num_shifts, np.float32)
        for charge in range(1, num_shifts):
            mass_diff[charge] = precursor_mass_diff / charge

        # Find the matching peaks between both spectra.
        peak_match_scores, peak_match_idx = [], []
        for peak_index, (peak_mz, peak_intensity) in enumerate(
                zip(spec.mz, spec.intensity)
        ):
            # Advance while there is an excessive mass difference.
            for cpi in range(num_shifts):
                while other_peak_index[cpi] < len(spec_other.mz) - 1 and (
                        peak_mz - fragment_mz_tolerance
                        > spec_other.mz[other_peak_index[cpi]] + mass_diff[cpi]
                ):
                    other_peak_index[cpi] += 1
            # Match the peaks within the fragment mass window if possible.
            for cpi in range(num_shifts):
                index = 0
                other_peak_i = other_peak_index[cpi] + index
                while (
                        other_peak_i < len(spec_other.mz)
                        and abs(
                    peak_mz - (spec_other.mz[other_peak_i] + mass_diff[cpi])
                )
                        <= fragment_mz_tolerance
                ):
                    peak_match_scores.append(
                        peak_intensity * spec_other.intensity[other_peak_i]
                    )
                    peak_match_idx.append((peak_index, other_peak_i))
                    index += 1
                    other_peak_i = other_peak_index[cpi] + index

        score, peak_matches = 0.0, []
        if len(peak_match_scores) > 0:
            # Use the most prominent peak matches to compute the score (sort in
            # descending order).
            peak_match_scores_arr = np.asarray(peak_match_scores)
            peak_match_order = np.argsort(peak_match_scores_arr)[::-1]
            peak_match_scores_arr = peak_match_scores_arr[peak_match_order]
            peak_match_idx_arr = np.asarray(peak_match_idx)[peak_match_order]
            peaks_used, other_peaks_used = set(), set()
            for peak_match_score, peak_i, other_peak_i in zip(
                    peak_match_scores_arr,
                    peak_match_idx_arr[:, 0],
                    peak_match_idx_arr[:, 1],
            ):
                if (
                        peak_i not in peaks_used
                        and other_peak_i not in other_peaks_used
                ):
                    score += peak_match_score
                    # Save the matched peaks.
                    peak_matches.append((peak_i, other_peak_i))
                    # Make sure these peaks are not used anymore.
                    peaks_used.add(peak_i)
                    other_peaks_used.add(other_peak_i)

        # return cosine similarity score and peak_match indices
        return score, peak_matches

    # create undirected connected graph for topological sort using mst edges
    G = nx.Graph()

    # add edges
    for edge in mst.edges(data=True):
        spec1, spec2, weight = edge
        G.add_edge(int(spec1.get('scans')), int(spec2.get('scans')))

    bfs = nx.bfs_tree(G, list(G.nodes)[0])  # change mst to bfs to set a starting node/change to acylic graph

    # topological sort
    topo_sorted_nodes = list(nx.topological_sort(bfs))

    print(topo_sorted_nodes)

    # store the alignment results between connected nodes in topological sorting in a dictionary
    alignment_results = {}

    idx = 0
    while (idx < len(topo_sorted_nodes) - 1):
        # get indices
        spec1_idx = topo_sorted_nodes[idx]
        spec2_idx = topo_sorted_nodes[idx + 1]
        spec1 = spec_dic[spec1_idx]
        spec2 = spec_dic[spec2_idx]
        # calculate cosine scores and peak matches
        score, peak_matches = _cosine_fast(spec1, spec2, 0.1, True)
        # store results
        alignment_results[spec1_idx, spec2_idx] = peak_matches
        print(f"scans: {spec1_idx} - {spec2_idx} and match indices: {peak_matches}")
        idx = idx + 1

    # xianghu's function needed for get_sets function
    def peak_tuple_to_dic(peakmatches):
        dic = {}
        for peakmatch in peakmatches:
            dic[peakmatch[0]] = peakmatch[1]
        return dic

    # peak tuple for the each peak in a set to store the scan number and peak idx
    PeakTuple = collections.namedtuple("PeakTuple", ["scan_num", "peak_idx"])

    # function to get sets given a path
    def get_sets(path):
        # initialize list of sets/lists
        final_match_list = []

        # create initial sets from the first alignment
        for match in alignment_results[path[0], path[1]]:
            tuple_1 = PeakTuple(path[0], match[0])
            tuple_2 = PeakTuple(path[1], match[1])
            final_match_list.append([tuple_1, tuple_2])

        idx = 1
        
        # iterate through the rest of the path
        while idx < len(path) - 1:
            peak_matches = alignment_results[path[idx], path[idx + 1]]
            peak_dic = peak_tuple_to_dic(peak_matches)
            temp_match_list = []

            # iterate through each set (stored as a list) in the list of sets
            for current_set in final_match_list:
                peak_to_find = list(current_set)[-1].peak_idx
                peak_scan = list(current_set)[-1].scan_num

                if (peak_to_find in peak_dic and peak_scan == path[idx]):
                    value = peak_dic[peak_to_find]
                    current_set.append(PeakTuple(path[idx + 1], value))
                    # keep track of peak pairs that were merged into an existing set
                    temp_match_list.append([PeakTuple(path[idx], peak_to_find), PeakTuple(path[idx + 1], value)])

            # append pairs that were not added to existing sets
            for peak in peak_matches:
                tuple_1 = PeakTuple(path[idx], peak[0])
                tuple_2 = PeakTuple(path[idx + 1], peak[1])
                new_pair = [tuple_1, tuple_2]
                if new_pair not in temp_match_list:
                    final_match_list.append(new_pair)

            idx += 1

        return final_match_list

    # get sets
    path = topo_sorted_nodes
    sets = get_sets(path)

    # function to get pairs that are not directly connected by an edge
    def new_matches (topo_sorted_nodes):
        idx1 = 0
        final_list_tuples = []

        while (idx1 < len(topo_sorted_nodes) - 2):
            # get indices
            spec1_idx = topo_sorted_nodes[idx1]
            idx2 = idx1 + 2 # +2 to skip over the directly connected node
            
            while (idx2 < len(topo_sorted_nodes)):
                spec2_idx = topo_sorted_nodes[idx2] 
                spec1 = spec_dic[spec1_idx]
                spec2 = spec_dic[spec2_idx]
                # calculate cosine scores and peak matches
                score, peak_matches = _cosine_fast(spec1, spec2, 0.1, True)
                # store results
                # alignment_results[spec1_idx, spec2_idx] = peak_matches
                # store results in a list of list of peak tuples
                for match in peak_matches:
                    tuple_1 = PeakTuple(spec1_idx, match[0])
                    # print(match[0])
                    tuple_2 = PeakTuple(spec2_idx, match[1])
                    final_list_tuples.append([tuple_1, tuple_2])

                print(f"scans: {spec1_idx} - {spec2_idx} and match indices: {peak_matches}")
                for match in peak_matches:
                    print(f"spec {spec1_idx} m/z: {spec1.mz[match[0]]}, spec {spec2_idx} m/z: {spec2.mz[match[1]]}")

                idx2 = idx2 + 1

            idx1 = idx1 + 1

        return final_list_tuples


    new_pairs_list = new_matches(topo_sorted_nodes)

    # function to merge pairs into transitively aligned sets or add them to the list of sets
    def add_pairs(transitive_sets, new_sets):
        peak_to_set = {}
    
        # Initialize peak_to_set dictionary to map each peak to its corresponding set index
        for set_idx, current_set in enumerate(transitive_sets):
            for peak in current_set:
                peak_to_set[(peak.scan_num, peak.peak_idx)] = set_idx

        for pair in new_sets:
            peak1 = (pair[0].scan_num, pair[0].peak_idx)
            peak2 = (pair[1].scan_num, pair[1].peak_idx)
            
            set1 = peak_to_set.get(peak1)
            set2 = peak_to_set.get(peak2)
            
            if set1 is not None and set2 is None:
                if pair[1].scan_num not in [peak.scan_num for peak in transitive_sets[set1]]:
                    transitive_sets[set1].append(pair[1])
                    peak_to_set[peak2] = set1
            
            elif set1 is None and set2 is None:
                new_set_idx = len(transitive_sets)
                transitive_sets.append([pair[0], pair[1]])
                peak_to_set[peak1] = new_set_idx
                peak_to_set[peak2] = new_set_idx

        # remove empty sets and duplicates
        transitive_sets = [list({(peak.scan_num, peak.peak_idx) for peak in t_set}) for t_set in transitive_sets if t_set]

        # convert sets back to list of PeakTuples
        transitive_sets = [[PeakTuple(scan_num, peak_idx) for scan_num, peak_idx in t_set] for t_set in transitive_sets]

        return transitive_sets

    final_sets = add_pairs(sets, new_pairs_list)

    # reorder spec_dic to be in topological order
    spec_dic = {key: spec_dic[key] for key in path}

    new_spec_dic = {}
    # convert spectrum tuple to dictionary
    for key, s in spec_dic.items():
        new_spec_dic[int(s.scan)] = {
            "scan": int(s.scan),
            "precursor_mz": float(s.precursor_mz),
            "precursor_charge": int(s.precursor_charge),
            "mz": s.mz.tolist(),  # convert to list
            "intensity": s.intensity.tolist()  # convert to list
        }

    # convert sets (list of lists with named tuple) to list of tuples
    # function to convert peak tuple to tuple
    def peak_tuple_to_tuple(peak_tuple):
        return (int(peak_tuple.scan_num), int(peak_tuple.peak_idx))

    # change sets to list of lists, each element in the list is a tuple
    json_sets = []
    for match_list in final_sets:
        json_match_list = [peak_tuple_to_tuple(peak_tuple) for peak_tuple in match_list]
        json_sets.append(json_match_list)

    return component, json_sets, new_spec_dic

# dash app to take input component number and process data
app = dash.Dash(__name__)

SAVE_PATH = os.path.join(os.getcwd(), 'saved.json')

@app.callback(
    Output('output-path', 'children'),
    [Input('process-button', 'n_clicks')],
    [State('task-id-input', 'value'), State('component-input', 'value')]
)
def process_and_save_json(n_clicks, task_id, component_number):
    if n_clicks is None or not task_id or not component_number:
        return dash.no_update
    
    preprocessed_spectra = get_data(int(component_number), task_id)
    component, json_sets, new_spec_dic = get_final_sets(int(component_number), preprocessed_spectra)

    # result_json = json.dumps({"component": component, "sets": json_sets, "spectra": new_spec_dic})

    # file_path = f'results_{task_id}_{component_number}.json'
    # with open(file_path, 'w') as f:
    #     f.write(result_json)

    # return f'JSON file saved at: {file_path}'

    with open(SAVE_PATH, 'w') as f:
        json.dump((component, json_sets, new_spec_dic), f)
    
    return f'JSON file saved at: {SAVE_PATH}'

app.layout = html.Div([
    dcc.Input(id='task-id-input', type='text', placeholder='Enter task ID'),
    dcc.Input(id='component-input', type='text', placeholder='Enter component number'),
    html.Button('Process and Save JSON', id='process-button'),
    html.Div(id='output-path')
])
    
if __name__ == '__main__':
    app.run(debug=True)