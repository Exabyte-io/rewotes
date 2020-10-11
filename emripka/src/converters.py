import pandas as pd
import json   

def csv_to_json(MP_output_dir, MP_json_dir, fname):
    """ Simple function to convert data from https://materialsproject.org/ to JSON format.

    Arguments:
        MP_output_dir (string):
            Directory where the file to convert exists.
                Example: "training/materialsproject_output/"

        MP_json_dir (string):
            Directory where the JSON file will be created.
                Example: "training/materialsproject_json/"

        fname (string):
            This is the file you wish to convert, without the file end.
                Example: "Si"

    Returns:
        Creates a new .json file in the MP_json_dir from the .csv file in the
        MP_output_dir.

    """
    # to-do: default directories for no input

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
        "formation_energy__eV",
        "E_above_hull__eV",
        "band_gap__eV",
        "has_bandstructure",
        "volume",
        "Nsites",
        "theoretical",
        "count",
        "density__gm_per_cc",
        "crystal_system",
    ]
    data_df = data_df.set_index('material_ID')

    # convert dataframe to dictionary of dictionaries
    data_dic = data_df.to_dict('index')

    # dump dictionary to a json file
    with open(f"{MP_json_dir}/{fname}.json","w") as file:
        json.dump(data_dic,file,indent=4)

def create_non_numeric_map(data_dict, non_numeric_key):
    unique_keys_map = dict()
    ii = 0
    for value in data_dict[non_numeric_key]:
        if value not in list(unique_keys_map.keys()):
            unique_keys_map[value] = ii 
            ii += 1
    return unique_keys_map
