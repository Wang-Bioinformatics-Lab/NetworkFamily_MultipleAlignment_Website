import networkx as nx
import collections
import numpy as np
import math
from typing import List, Tuple, Dict, Set
from itertools import combinations
import requests
import io
import pandas as pd
import argparse
import json

SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["scan", "precursor_mz", "precursor_charge", "mz", "intensity"]
)

# peak tuple for the each peak in a set to store the scan number and peak idx
PeakTuple = collections.namedtuple("PeakTuple", ["scan_num", "peak_idx"])

# xianghu's functions
def norm_intensity(intensity):
    return np.copy(intensity) / np.linalg.norm(intensity)

def read_mgf(file_data):
    spectra = []  # List to store all spectra
    spectrum = None  # Current spectrum being processed
    
    for line in file_data:
        line = line.strip()
        
        if not line:
            # Skip empty lines
            continue
        
        if line == 'BEGIN IONS':
            spectrum = {'peaks': [], 'm/z array': [], 'intensity array': []}  # Initialize a new spectrum
        elif line == 'END IONS':
            spectra.append(spectrum)  # Add the completed spectrum to the list
            spectrum = None  # Reset for the next spectrum
        elif spectrum is not None:
            if '=' in line:  # Property line
                key, value = line.split('=', 1)  # Split at the first '=' only
                spectrum[key.lower()] = value  # Convert key to lowercase for consistency
            else:  # Peak line
                try:
                    mz, intensity = line.split()
                    spectrum['peaks'].append((float(mz), float(intensity)))
                    spectrum['m/z array'].append(float(mz))
                    spectrum['intensity array'].append(float(intensity))
                except ValueError:
                    # Optionally, log a warning about the malformed line
                    print(f"Warning: Skipping malformed line: '{line}'")
                    continue
    
    return spectra

def mgf_processing(spectra):
    spec_dic = {}
    for spectrum in spectra:
        mz_array = spectrum['m/z array']
        intensity_array = spectrum['intensity array']
        filtered_mz = []
        filtered_intensities = []
        precursor_value = float(spectrum['pepmass'])
        charge = int(spectrum['charge'].rstrip('+'))
        scans = int(spectrum['scans'])
        for i, mz in enumerate(mz_array):
            peak_range = [j for j in range(len(mz_array)) if abs(mz_array[j] - mz) <= 25]
            sorted_range = sorted(peak_range, key=lambda j: intensity_array[j], reverse=True)
            if i in sorted_range[:6]:
                if abs(mz - precursor_value) > 17:
                    filtered_mz.append(mz)
                    filtered_intensities.append(intensity_array[i])
        filtered_intensities = [math.sqrt(x) for x in filtered_intensities]
        spec_dic[scans] = SpectrumTuple(scans, precursor_value, charge, filtered_mz,
                                                       norm_intensity(filtered_intensities))
    return spec_dic

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

# xianghu's function needed for get_sets function
def peak_tuple_to_dic(peakmatches):
    dic = {}
    for peakmatch in peakmatches:
        dic[peakmatch[0]] = peakmatch[1]
    return dic


# function to get data from csv and mgf files
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
    # spectrums = list(load_from_mgf(mgf_data))
    spectra = read_mgf(mgf_data)

    # preprocess spectra
    spec_dic = mgf_processing(spectra)
    
    # process DataFrame
    df_comp = df[df['ComponentIndex'] == component]
    # print(df_comp)
    df_comp = df_comp[['CLUSTERID1', 'CLUSTERID2', 'ComponentIndex', 'Cosine']]
    scan_numbers1 = set(df_comp['CLUSTERID1'].astype(int).tolist())
    scan_numbers2 = set(df_comp['CLUSTERID2'].astype(int).tolist())
    scan_numbers = set(scan_numbers2) | set(scan_numbers1)

    filtered_spec_dic ={}

    # get spectra that match the scan numbers
    for scan in scan_numbers:
        if scan in spec_dic:
            filtered_spec_dic[scan] = spec_dic[scan] 
    
    return filtered_spec_dic, df_comp, scan_numbers
    
# function to get mst and do topo sort
def get_topo_path(filtered_spec_dic, df_comp, scan_nums):
    pairwise_cosine_scores = []
    
    # Get all pairs of scan numbers
    scan_pairs = combinations(scan_nums, 2)
    
    for scan1, scan2 in scan_pairs:
        spec1 = filtered_spec_dic[scan1]
        spec2 = filtered_spec_dic[scan2]

        # Calculate cosine score
        cosine_score, pairs = _cosine_fast(spec1, spec2, 0.1, True)

        # Store the results as a tuple (scan1, scan2, cosine_score)
        pairwise_cosine_scores.append((scan1, scan2, cosine_score))


    # make an undirected graph for mst
    graph = nx.Graph()

    # add nodes
    for i, s in filtered_spec_dic.items():
        graph.add_node(i, scan = s.scan)

    # print(scores)

    # add edges/modified cosine similarities
    for (spec1, spec2, score) in pairwise_cosine_scores:
        if spec1 != spec2:
            graph.add_edge(spec1, spec2, weight = 1 - score)
            
    # kruskal's min spanning tree algorithm
    mst = nx.minimum_spanning_tree(graph)

    # output mst edges
    # print("MST Edges:")
    # for edge in mst.edges(data=True):
    #     spec1, spec2, weight = edge
    #     print(f"(Scan {spec1}) - "
    #         f"(Scan {spec2}), "
    #         f"Similarity Score: {1 - weight['weight']}")
    #         # f"Similarity Score: {weight['Cosine']}")

    # create undirected connected graph for topological sort using mst edges
    G = nx.Graph()

    # add edges
    for edge in mst.edges(data=True):
        spec1, spec2, data = edge
        G.add_edge(spec1, spec2)

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
        spec1 = filtered_spec_dic[spec1_idx]
        spec2 = filtered_spec_dic[spec2_idx]
        # calculate cosine scores and peak matches
        score, peak_matches = _cosine_fast(spec1, spec2, 0.1, True)
        # store results
        alignment_results[spec1_idx, spec2_idx] = peak_matches
        print(f"scans: {spec1_idx} - {spec2_idx} and match indices: {peak_matches}")
        idx = idx + 1

    return topo_sorted_nodes, alignment_results

# function to get sets given a path
def get_sets(path, alignment_results):
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

# function to get pairs that are not directly connected by an edge
def new_matches (topo_sorted_nodes, filtered_spec_dic):
    idx1 = 0
    final_list_tuples = []

    while (idx1 < len(topo_sorted_nodes) - 2):
        # get indices
        spec1_idx = topo_sorted_nodes[idx1]
        idx2 = idx1 + 2 # +2 to skip over the directly connected node
        
        while (idx2 < len(topo_sorted_nodes)):
            spec2_idx = topo_sorted_nodes[idx2] 
            spec1 = filtered_spec_dic[spec1_idx]
            spec2 = filtered_spec_dic[spec2_idx]
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

            # print(f"scans: {spec1_idx} - {spec2_idx} and match indices: {peak_matches}")
            # for match in peak_matches:
            #     print(f"spec {spec1_idx} m/z: {spec1.mz[match[0]]}, spec {spec2_idx} m/z: {spec2.mz[match[1]]}")

            idx2 = idx2 + 1

        idx1 = idx1 + 1

    return final_list_tuples

# function to merge pairs into transitively aligned sets or add them to the list of sets
def add_pairs(transitive_sets, new_sets):
    peak_to_set = {}

    # Initialize peak_to_set dictionary to map each peak to its corresponding set index
    for set_idx, current_set in enumerate(transitive_sets):
        for peak in current_set:
            peak_to_set[(peak.scan_num, peak.peak_idx)] = set_idx

    sets_to_remove = set()

    for pair in new_sets:
        peak1 = (pair[0].scan_num, pair[0].peak_idx)
        peak2 = (pair[1].scan_num, pair[1].peak_idx)
        
        set1 = peak_to_set.get(peak1)
        set2 = peak_to_set.get(peak2)
        
        # if the first peak is in a set but the second is not
        if set1 is not None and set2 is None:
            # add the peak if the scan num is not already in the set
            if pair[1].scan_num not in [peak.scan_num for peak in transitive_sets[set1]]:
                transitive_sets[set1].append(pair[1])
                peak_to_set[peak2] = set1

        # if the second peak is in a set but the first is not
        elif set1 is None and set2 is not None:
            # add the peak if the scan num is not already in the set
            if pair[0].scan_num not in [peak.scan_num for peak in transitive_sets[set2]]:
                transitive_sets[set2].append(pair[0])
                peak_to_set[peak1] = set2
        
        elif set1 is None and set2 is None:
            new_set_idx = len(transitive_sets)
            transitive_sets.append([pair[0], pair[1]])
            peak_to_set[peak1] = new_set_idx
            peak_to_set[peak2] = new_set_idx
        
        # if both peaks are in sets already
        elif set1 is not None and set2 is not None:
            # check if any of the scan nums in set1 are in set 2 and if they are not, merge the sets
            set1_scans = {peak.scan_num for peak in transitive_sets[set1]}
            set2_scans = {peak.scan_num for peak in transitive_sets[set2]}
            if not set1_scans.intersection(set2_scans):
                transitive_sets[set1].extend(transitive_sets[set2])
                # update peak_to_set mapping for the merged set
                for peak in transitive_sets[set2]:
                    peak_to_set[(peak.scan_num, peak.peak_idx)] = set1
        
                # remove set2
                sets_to_remove.add(set2)

    # get final sets
    transitive_sets = [s for i, s in enumerate(transitive_sets) if i not in sets_to_remove]

    # remove empty sets and duplicates
    transitive_sets = [list({(peak.scan_num, peak.peak_idx) for peak in t_set}) for t_set in transitive_sets if t_set]

    # convert sets back to list of PeakTuples
    transitive_sets = [[PeakTuple(scan_num, peak_idx) for scan_num, peak_idx in t_set] for t_set in transitive_sets]

     # output sets
    # for i, s in enumerate(transitive_sets):
    #     print(f"set number: {i + 1}")
    #     for peak_tuple in s:
    #         mz = filtered_spec_dic[peak_tuple.scan_num].mz[peak_tuple.peak_idx]
    #         print(peak_tuple)
    #         print(mz)
    return transitive_sets

# function to convert peak tuple to tuple
def peak_tuple_to_tuple(peak_tuple):
    return (int(peak_tuple.scan_num), int(peak_tuple.peak_idx))

def data_for_json (final_sets, filtered_spec_dic, path):
    # reorder spec_dic to be in topological order
    filtered_spec_dic = {key: filtered_spec_dic[key] for key in path}

    new_spec_dic = {}

    # convert spectrum tuple to dictionary
    for key, s in filtered_spec_dic.items():
        new_spec_dic[int(s.scan)] = {
            "scan": int(s.scan),
            "precursor_mz": float(s.precursor_mz),
            "precursor_charge": int(s.precursor_charge),
            "mz": s.mz,
            "intensity": s.intensity.tolist() # convert to list
        }

    # change sets to list of lists, each element in the list is a tuple
    json_sets = []
    for match_list in final_sets:
        # convert sets (list of lists with named tuple) to list of tuples
        json_match_list = [peak_tuple_to_tuple(peak_tuple) for peak_tuple in match_list]
        json_sets.append(json_match_list)

    return json_sets, new_spec_dic

# def data_for_json(final_sets, filtered_spec_dic, path):
#     # Reorder spec_dic to be in topological order
#     filtered_spec_dic = {key: filtered_spec_dic[key] for key in path}

#     new_spec_dic = {}

#     # Create a dictionary to map peaks to sets
#     peak_to_set = {}
#     for set_idx, match_list in enumerate(final_sets):
#         for peak_tuple in match_list:
#             peak_to_set[peak_tuple] = set_idx

#     # Convert spectrum tuple to dictionary
#     for key, s in filtered_spec_dic.items():
#         # Get the top 10 peaks based on intensity
#         top_10_indices = np.argsort(s.intensity)[-10:][::-1]
#         top_10_mz = [s.mz[i] for i in top_10_indices]
#         top_10_intensity = [s.intensity[i] for i in top_10_indices]

#         # Add a tag for peaks that are part of a set
#         peaks_tagged = []
#         for i, (mz, intensity) in enumerate(zip(s.mz, s.intensity)):
#             tag = peak_to_set.get(PeakTuple(key, i), None)
#             peaks_tagged.append({
#                 "mz": mz,
#                 "intensity": intensity,
#                 "in_set": tag is not None,
#                 "set_id": tag if tag is not None else -1  # -1 indicates no set
#             })

#         new_spec_dic[int(s.scan)] = {
#             "scan": int(s.scan),
#             "precursor_mz": float(s.precursor_mz),
#             "precursor_charge": int(s.precursor_charge),
#             "mz": s.mz,
#             "intensity": s.intensity.tolist(), # convert to list
#             "top_10_mz": top_10_mz,
#             "top_10_intensity": top_10_intensity,
#             "peaks_tagged": peaks_tagged,
#         }

#     # Change sets to list of lists, each element in the list is a tuple
#     json_sets = []
#     for match_list in final_sets:
#         # Convert sets (list of lists with named tuple) to list of tuples
#         json_match_list = [peak_tuple_to_tuple(peak_tuple) for peak_tuple in match_list]
#         json_sets.append(json_match_list)

#     return json_sets, new_spec_dic


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some spectrum data and generate transitively aligned sets.')
    
    # input for get_data
    parser.add_argument('-component', type=int, required=True, help='Component index to filter spectra')
    parser.add_argument('-task_id', type=str, required=True, help='Task ID for retrieving data from GNPS')
    # output file name
    parser.add_argument('-output', type=str, help='Output file name')

    args = parser.parse_args()

    # get data using the task_id and component
    filtered_spec_dic, df_comp, scan_numbers = get_data(args.component, args.task_id)
    
    # get topological path and alignment results
    topo_sorted_nodes, alignment_results = get_topo_path(filtered_spec_dic, df_comp, scan_numbers)
    
    # get sets using the topological path
    transitive_sets = get_sets(topo_sorted_nodes, alignment_results)
    
    # Generate new matches for nodes not directly connected by edges
    new_pairs = new_matches(topo_sorted_nodes, filtered_spec_dic)
    
    # Merge new pairs into existing sets
    final_sets = add_pairs(transitive_sets, new_pairs)

    json_sets, new_spec_dic = data_for_json(final_sets, filtered_spec_dic, topo_sorted_nodes)
    

    output_path = args.output if args.output.endswith('.json') else f"{args.output}.json"
    with open(output_path, 'w') as f:
        json.dump((json_sets, new_spec_dic), f)
