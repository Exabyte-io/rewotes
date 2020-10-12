import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import json
from os import listdir
from os.path import isfile, join
from . import converters as converters
from . import stoichiometry as stoichiometry

import os
this_dir, this_filename = os.path.split(__file__)

def create_id():
    """
    Creates a unique id for each user input TrainingData object.

    Returns:
        (str): of format "user-#####"

    To-do: need to make this a unique value wrt all user data made and stored
    """
    return "user-%05d" % np.random.randint(0,99999)  

class TrainingData:
    """
    User input training data class. A formula and band_gap are required, with all other input parameters
    optional, but will strengthen the model if using just user training data to predict.

    Arguments:
        formula (str)
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
    def __init__(self,formula,band_gap,spacegroup="",formation_energy=0.0,E_above_hull=0.0,has_bandstructure=False,
                    volume=0.0,Nsites=0,theoretical=False,count=0.0,density=0.0,crystal_system=""):
        self.ID = create_id() 
        self.formula = formula
        self.spacegroup = spacegroup 
        self.formation_energy = formation_energy 
        self.E_above_hull = E_above_hull 
        self.band_gap = band_gap 
        self.has_bandstructure = has_bandstructure 
        self.volume = volume 
        self.Nsites = Nsites 
        self.theoretical = theoretical 
        self.count = count 
        self.density = density 
        self.crystal_system = crystal_system 
        self.data_dict = dict() 

    def make_data_dict(self):
        """
        Populate's the object's data_dict a dictionary from input parameters 
        to match the structure found in materialsproject.org's json files.
        """
        data_dict = {
            self.ID: {
                "formula": self.formula,
                "spacegroup": self.spacegroup,
                "formation_energy__eV": self.formation_energy,
                "E_above_hull__eV": self.E_above_hull,
                "band_gap__eV": self.band_gap,
                "has_bandstructure": self.has_bandstructure,
                "volume": self.volume,
                "Nsites": self.Nsites,
                "theoretical": self.theoretical,
                "count": self.count,
                "density__gm_per_cc": self.density,
                "crystal_system": self.crystal_system 
            }
        }
        self.data_dict = data_dict

    def show_data(self):
        """
        Helper function for user to see the parameters they have input. 
        """
        self.make_data_dict()
        print(self.data_dict)

    def store_data(self):
        """ 
        Stores the user training data to the data directory in the 
        user_data.json file for use in training the model.
        """
        self.make_data_dict()

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
        periodic_table_obj (PeriodicTable object)
        use_database_data (bool): decides useage of the package
            True: the database data is used to train the model, along with any input TrainingData objects 
            False: the database data is not used to train the model, and only input TrainingData is
    To-do:
        - Move the user_data.json file to /data/training/user_data_json/ and rename to the date and time
          of the prediciton after predictions
    """
    def __init__(self, periodic_table_obj, use_database_data):
        self.use_database_data = use_database_data
        self.csv_path = this_dir+"/data/training/materialsproject_output/"
        self.json_path = this_dir+"/data/training/materialsproject_json/"
        self.periodic_table_obj = periodic_table_obj
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
                converters.csv_to_json(self.csv_path,self.json_path,fname=csv_file)
                print(f"created {csv_file}.json")

    def get_stored_data(self):
        """
        Accesses the stored json file data and populates it in self.data_dict.
            - sets the molecular_weight parameter as well
        """
        if self.use_database_data:
            json_files = [f for f in listdir(self.json_path) if isfile(join(self.json_path, f))]
        else:
            json_files = [self.json_path+"user_data.json"]
        training_compounds = [f.split(".json")[0] for f in json_files]

        for training_compound in training_compounds:
            with open(f"{self.json_path}/{training_compound}.json") as fname:
                training_compound_results = dict(json.load(fname))
                for ID, training_compound_result in training_compound_results.items():
                    new_result = training_compound_result
                    new_result["molecular_weight"] = stoichiometry.get_molecular_weight(new_result["formula"], self.periodic_table_obj)  
                    self.data_dict[ID] = new_result

        self.data_IDs = list(self.data_dict.keys())

    def set_stoichiometry(self):
        """
        Sets the stoichiometry parameters for the material.
            Example: { "Cu": 0.4, "S": 0.6 }
        """
        for ID, result in self.data_dict.items():
            formula = result["formula"]
            norm_stoichiometry = stoichiometry.get_norm_stoichiomertry(formula)  
            self.data_dict[ID]["stoichiometry"] = norm_stoichiometry 

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
        self.non_element_keys = ["band_gap__eV",]
        for param in material_training_params:
            self.non_element_keys.append(param)
        self.non_element_keys.append("molecular_weight")
        self.populate_data_dict_clean()

        # create dataframe from data_dict_clean 
        self.dataframe = pd.DataFrame(self.data_dict_clean) 

        # dropping rows which have a bandgap of zero; not sure if I want to do this...
        self.dataframe = self.dataframe[self.dataframe['band_gap__eV'] != 0]

    def populate_data_dict_clean(self):
        """
        Populates data_dict_clean with training params only in self.non_element_keys.
        """
        for non_element_key in self.non_element_keys:
            self.data_dict_clean[non_element_key] = [result[non_element_key] for (ID, result) in self.data_dict.items() ] 

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

    def get_train_test_splits(self):
        """
        Helper function used to extract the correctly formatted data for use in training
        the model.

        Returns:
            X_train (arr)
            X_test (arr)
            y_train (arr)
            y_test (arr)
        """
        X_keys = list(self.dataframe.keys())[2:]

        X = np.asarray(self.dataframe[X_keys])
        y = np.asarray(self.dataframe['band_gap__eV'])
                
        return train_test_split(X, y, test_size=0.33, shuffle= True) 
