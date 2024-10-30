import collections
import math
import requests

SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["scan", "precursor_mz", "precursor_charge", "mz", "intensity"]
)

def get_usi_data(usi_url):
    response = requests.get(usi_url)
    response.raise_for_status()
    return response.json()

def usi_processing(usi_string):
    spec_dic = {}
    usi_list = usi_string.split(',')
    
    for usi_url in usi_list:
        usi_data = get_usi_data(usi_url.strip())

        precursor_mz = float(usi_data.get("precursor_mz", 0.0))
        precursor_charge = int(usi_data.get("precursor_charge", 0))
        scan = int(usi_data.get("scan", 0))
        peaks = usi_data.get("peaks", [])

        mz_array = [peak[0] for peak in peaks]
        intensity_array = [peak[1] for peak in peaks]