import pandas as pd
import json   
import os
this_dir, this_filename = os.path.split(__file__)

def csv_to_json(fname):
    """ Simple function to convert data from https://materialsproject.org/ to JSON format.

    Args:
        fname (str)
            This is the file you wish to convert, without the file end.
                Example: "Si"

    Returns:
        Creates a new .json file in the MP_json_dir from the .csv file in the
        MP_output_dir.

    """
    MP_output_dir = this_dir+"/data/training/materialsproject_output/"
    MP_json_dir = this_dir+"/data/training/materialsproject_json/"

    # remove file end if still input
    if ".csv" in fname:
        fname = fname.split(".csv")[0]

    # use pandas to read csv into dataframe
    data_df = pd.read_csv(f"{MP_output_dir}/{fname}.csv")

    # clean up dataframe
    data_df = data_df.drop('Unnamed: 13', axis=1)
    data_df.columns = [
        "material_ID",
        "formula",
        "spacegroup",
        "formation_energy",
        "E_above_hull",
        "band_gap",
        "has_bandstructure",
        "volume",
        "Nsites",
        "theoretical",
        "count",
        "density",
        "crystal_system",
    ]
    data_df = data_df.set_index('material_ID')

    # convert dataframe to dictionary of dictionaries
    data_dic = data_df.to_dict('index')

    # dump dictionary to a json file
    with open(f"{MP_json_dir}/{fname}.json","w") as file:
        json.dump(data_dic,file,indent=4)

def create_non_numeric_map(data_dict, non_numeric_key):
    """
    Creates a dictionary used to map non-numeric training parameters to 
    unique identifier integer value.

    Args:
        data_dict (dict)
        non_numeric_key (str)

    Returns:
        unique_keys_map (dict)
    """
    unique_keys_map = dict()
    ii = 0
    for value in data_dict[non_numeric_key]:
        if value not in list(unique_keys_map.keys()):
            unique_keys_map[value] = ii 
            ii += 1
    return unique_keys_map
