import collections
import math
import requests
import urllib.parse

from alignment import norm_intensity

SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["scan", "precursor_mz", "precursor_charge", "mz", "intensity"]
)

def get_usi_url(usi):
    # base url
    base_url = "https://metabolomics-usi.gnps2.org/json/?usi1="

    # encode input usi
    encoded_usi = urllib.parse.quote(usi)
    # full url
    usi_url = f"{base_url}{encoded_usi}"

    response = requests.get(usi_url)
    response.raise_for_status()

    return response.json()

def parse_usi(usi):
    try:
        task_id = usi.split("TASK-")[1].split("-")[0]
    except IndexError:
        task_id = "unknown_task"

    try:
        scan = int(usi.split("scan:")[1].split("&")[0])
    except (IndexError, ValueError):
        scan = None

    return task_id, scan

def usi_processing(usi_string):
    spec_dic = {}
    # usi_list = usi_string.split(',')
    usi_list = usi_string.splitlines()
    scan_names = []
    
    for usi in usi_list:
        usi_data = get_usi_url(usi.strip())

        precursor_mz = float(usi_data.get("precursor_mz", 0.0))
        precursor_charge = int(usi_data.get("precursor_charge", 0))
        task_id, scan = parse_usi(usi)  # get task ID and scan number to create unique scan

        # for unique taskid-scan names
        # scan_name = f"{task_id}_scan{scan}"
        # scan_names.append(scan_name)  # add scan names to the list

        # using scan numbers only for spectrum in one component
        scan_names.append(int(scan))

        peaks = usi_data.get("peaks", [])
        mz_array = [peak[0] for peak in peaks]
        intensity_array = [peak[1] for peak in peaks]

        # from xianghu's mgf_processing function
        filtered_mz = []
        filtered_intensities = []
        for i, mz in enumerate(mz_array):
            peak_range = [j for j in range(len(mz_array)) if abs(mz_array[j] - mz) <= 25]
            sorted_range = sorted(peak_range, key=lambda j: intensity_array[j], reverse=True)
            
            if i in sorted_range[:6] and abs(mz - precursor_mz) > 17:
                filtered_mz.append(mz)
                filtered_intensities.append(intensity_array[i])

        filtered_intensities = [math.sqrt(x) for x in filtered_intensities]
        # spec_dic[scan] = SpectrumTuple(scan_name, precursor_mz, precursor_charge, filtered_mz, 
        #                                norm_intensity(filtered_intensities))

        spec_dic[scan] = SpectrumTuple(int(scan), precursor_mz, precursor_charge, filtered_mz, 
                                       norm_intensity(filtered_intensities))
        
    return spec_dic, scan_names

def get_fbmn_json (task_id, cluster_num):
    url = f"https://metabolomics-usi.gnps2.org/json/?usi1=mzspec:GNPS2:TASK-{task_id}-nf_output/clustering/spectra_reformatted.mgf:scan:{cluster_num}"

    response = requests.get(url)
    response.raise_for_status()

    return response.json()

def process_fbmn_input (task_id, cluster_string):
    spec_dic = {}
    cluster_list = cluster_string.split(',')
    scan_names = []
    
    for cluster in cluster_list:
        cluster_data = get_fbmn_json(task_id, cluster.strip())

        precursor_mz = float(cluster_data.get("precursor_mz", 0.0))
        precursor_charge = int(cluster_data.get("precursor_charge", 0))

        # using scan numbers only for spectrum in one component
        scan_names.append(int(cluster))

        peaks = cluster_data.get("peaks", [])
        mz_array = [peak[0] for peak in peaks]
        intensity_array = [peak[1] for peak in peaks]

        # from xianghu's mgf_processing function
        filtered_mz = []
        filtered_intensities = []
        for i, mz in enumerate(mz_array):
            peak_range = [j for j in range(len(mz_array)) if abs(mz_array[j] - mz) <= 25]
            sorted_range = sorted(peak_range, key=lambda j: intensity_array[j], reverse=True)
            
            if i in sorted_range[:6] and abs(mz - precursor_mz) > 17:
                filtered_mz.append(mz)
                filtered_intensities.append(intensity_array[i])

        filtered_intensities = [math.sqrt(x) for x in filtered_intensities]
    
        spec_dic[int(cluster)] = SpectrumTuple(int(cluster), precursor_mz, precursor_charge, filtered_mz, 
                                       norm_intensity(filtered_intensities))
        
    return spec_dic, scan_names

def get_scans (task_id, component):
    # get the json data
    url = f"https://gnps2.org/result?task={task_id}&viewname=clustersummary&json"
    response = requests.get(url)
    response.raise_for_status()
    json_data = response.json()

    # get corresponding cluster indices for the component
    # cluster_indices = []
    # for cluster in json_data.get("clusters", []):
    #     if cluster.get("componentindex") == component:
    #         cluster_indices.append(cluster.get("clusterindex"))

    cluster_indices = [cluster["cluster index"] for cluster in json_data if cluster.get("component") == str(component)]


    # change cluster indices into a string
    cluster_string = ",".join(map(str, cluster_indices))

    spec_dic, scan_names = process_fbmn_input(task_id, cluster_string)

    return spec_dic, scan_names