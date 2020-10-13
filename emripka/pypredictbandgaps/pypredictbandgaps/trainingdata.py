import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import json
from os import listdir
from os.path import isfile, join
from . import converters as converters

from datetime import datetime
import os
this_dir, this_filename = os.path.split(__file__)

import pymatgen as mg
from pymatgen.ext.matproj import MPRester

def create_id():
    """
    Creates a unique id for each user input TrainingData object.

    Returns:
        (str): of format "user_ID_YYYY_MM_DD_HH_MM_SS"
    """
    date_time_str = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    ID = "%05d" % np.random.randint(0,99999)  
    return f"user_{ID}_{date_time_str}"   

class TrainingData:
    """
    User input training data class. A formula and band_gap are required, with all other input parameters
    optional, but will strengthen the model if using just user training data to predict.

    Args:
        formula (str)
    
    Kwargs:
        band_gap (float)
        spacegroup (str)
        formation_energy (float)
        E_above_hull (float)
        has_bandstructure (bool)
        volume (float)
        Nsites (int)
        theoretical (bool)
        density (float)
        crystal_system (str)
    """
    def __init__(self, formula, **kwargs):
        self.ID = create_id() 
        self.formula = formula
        self.features = dict(**kwargs)
        self.features["formula"] = formula
        self.data_dict = { self.ID: self.features } 

    def store_data(self):
        """ 
        Stores the user training data to the data directory in the 
        user_data.json file for use in training the model.
        """
        json_path = this_dir+"/data/training/materialsproject_json/"
        # open the existing user_data.json file
        with open(f"{json_path}/user_data.json","r") as fname:
            try:
                user_data = dict(json.load(fname))
            except:
                user_data = dict() 
        fname.close()

        user_data[self.ID] = self.data_dict[self.ID]
        with open(f"{json_path}/user_data.json","w") as fname:
            json.dump(user_data,fname,indent=4)
        fname.close()

class BandGapDataset:
    """
    Class to create a house the training data for the model.

    Arguments:
        use_database_data (bool): decides useage of the package
            True: the database data is used to train the model, along with any input TrainingData objects 
            False: the database data is not used to train the model, and only input TrainingData is
    """
    def __init__(self, use_database_data):
        self.use_database_data = use_database_data
        self.csv_path = this_dir+"/data/training/materialsproject_output/"
        self.json_path = this_dir+"/data/training/materialsproject_json/"
        self.data_dict = dict()
        self.data_IDs = list() 
        if self.use_database_data:
            self.convert_stored_csvs() 
        self.get_stored_data()
        self.set_stoichiometry()

    def convert_stored_csvs(self):
        """
        Converts the stored csv data to json files if this hasn't been previously done.
        """
        csv_files = [f for f in listdir(self.csv_path) if isfile(join(self.csv_path, f))]
        json_files = [f for f in listdir(self.json_path) if isfile(join(self.json_path, f))]

        csv_files_stripped = [f.split(".csv")[0] for f in csv_files]
        json_files_stripped = [f.split(".json")[0] for f in json_files]

        for csv_file in csv_files_stripped:
            if csv_file not in json_files_stripped:
                converters.csv_to_json(csv_file)
                print(f"created {csv_file}.json")

    def get_stored_data(self):
        """
        Accesses the stored json file data and populates it in self.data_dict.
            - sets the molecular_weight parameter as well
        """
        if self.use_database_data:
            json_files = [f for f in listdir(self.json_path) if isfile(join(self.json_path, f))]
        else:
            json_files = ["user_data.json"]
        training_compounds = [f.split(".json")[0] for f in json_files]

        for training_compound in training_compounds:
            with open(f"{self.json_path}/{training_compound}.json") as fname:
                training_compound_results = dict(json.load(fname))
                for ID, training_compound_result in training_compound_results.items():
                    new_result = training_compound_result
                    composition = mg.Composition(new_result["formula"])
                    new_result["molecular_weight"] = composition.weight  
                    self.data_dict[ID] = new_result

        self.data_IDs = list(self.data_dict.keys())

    def set_stoichiometry(self):
        """
        Sets the stoichiometry parameters for the material.
            Example: { "Cu": 0.4, "S": 0.6 }
        """
        for ID, result in self.data_dict.items():
            composition = mg.Composition(result["formula"])
            self.data_dict[ID]["stoichiometry"] = { element.value: composition.get_atomic_fraction(element) for element in composition }

class BandGapDataFrame:
    """
    Class which converts the band_gap_dataframe object to the correct structure needed
    to train the model.

    Arguments:
        data_dict (dict of dicts)
        symbols (list of str): symbols of all elements (from PeriodicTable object's symbols variable)
            Example: ["H", "He", ..., "Uuo"]
        material_training_params (list of str): contains the parameters which will be used to train the
            model (save the molecular_weight and symbols, as these will always be training params)
            Example: ["density","crystal_system"]
    """
    def __init__(self, data_dict, symbols, material_training_params):
        self.data_dict = data_dict
        self.symbols = symbols
        self.crystal_system_map = dict()
        self.spacegroup_map = dict()

        # create new dictionary
        self.data_dict_clean = {
            "ID": [ ID for (ID, result) in self.data_dict.items() ],
        }

        # select training params 
        self.non_element_keys = ["band_gap",]
        for param in material_training_params:
            if type(param) == str:
                self.non_element_keys.append(param)
        self.non_element_keys.append("molecular_weight")
        self.populate_data_dict_clean()

        # create dataframe from data_dict_clean 
        self.dataframe = pd.DataFrame(self.data_dict_clean) 

        # dropping rows which have a bandgap of zero; not sure if I want to do this...
        self.dataframe = self.dataframe[self.dataframe['band_gap'] != 0]

    def populate_data_dict_clean(self):
        """
        Populates data_dict_clean with training params only in self.non_element_keys.
        """
        for non_element_key in self.non_element_keys:
            self.data_dict_clean[non_element_key] = [ result[non_element_key] for (ID, result) in self.data_dict.items() ] 

        for symbol in self.symbols:
            self.data_dict_clean[symbol] = list() 

        self.populate_stoichiometry()

        if "crystal_system" in self.non_element_keys: 
            self.crystal_system_map = converters.create_non_numeric_map(self.data_dict_clean, "crystal_system")
            self.data_dict_clean["crystal_system"] = [self.crystal_system_map[value] for value in self.data_dict_clean["crystal_system"] ] 

        if "spacegroup" in self.non_element_keys: 
            self.spacegroup_map = converters.create_non_numeric_map(self.data_dict_clean, "spacegroup")
            self.data_dict_clean["spacegroup"] = [self.spacegroup_map[value] for value in self.data_dict_clean["spacegroup"] ] 
    
    def populate_stoichiometry(self):
        """
        Sets the stoichiometry parameters for the material for each element in the 
        periodic table.
        """
        for ID, result in self.data_dict.items():
            elements = list(result["stoichiometry"].keys())   
            for symbol in self.symbols:
                value = result["stoichiometry"][symbol] if symbol in elements else 0
                self.data_dict_clean[symbol].append(value)

    def get_train_test_splits(self, test_size=0.25):
        """
        Helper function used to extract the correctly formatted data for use in training
        the model.

        Arguments:
            test_size (float): optional input of test_size

        Returns:
            X_train (arr)
            X_test (arr)
            y_train (arr)
            y_test (arr)
        """
        X_keys = list(self.dataframe.keys())[2:]

        X = np.asarray(self.dataframe[X_keys])
        y = np.asarray(self.dataframe['band_gap'])
                
        return train_test_split(X, y, test_size=test_size, shuffle= True) 
