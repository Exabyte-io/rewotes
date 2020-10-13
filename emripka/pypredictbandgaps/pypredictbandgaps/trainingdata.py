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
    def __init__(self, material):
        self.material = material
        self.material_elements = list({ element.value: material.composition.get_atomic_fraction(element) for element in material.composition }.keys())
        self.data_dict = dict()
        self.get_data()

    def get_data(self):
        composition = mg.Composition(self.material.formula)
        chemical_system = composition.chemical_system

        with MPRester() as mpr:
            print(f"... getting mp ids for the chemical system: {chemical_system} ...")
            training_ids = mpr.get_materials_ids(chemical_system)        
            num_ids = len(training_ids)

            for idx, training_id in enumerate(training_ids):
                print(f" {training_id}: {idx+1} of {num_ids} mp ids")

                print(f"  ... getting structure ...")
                structure = mpr.get_structure_by_material_id(training_id)
                print(f"  ... getting band structure ...")
                try:
                    band_structure = mpr.get_bandstructure_by_material_id(training_id)
                except:
                    print("       failed")
                    continue

                tmp_data_dict = dict()
                if band_structure is not None: 
                    structure_dict = structure.as_dict()
                    band_structure_dict = band_structure.get_band_gap()
                    tmp_data_dict["formula"] = structure.composition.alphabetical_formula
                    for lattice_param in ["a","b","c","alpha","beta","gamma","volume"]:
                        tmp_data_dict[lattice_param] = structure_dict["lattice"][lattice_param]  
                    tmp_data_dict["band_gap"] = band_structure_dict["energy"]

                    composition = mg.Composition(tmp_data_dict["formula"])
                    tmp_data_dict["molecular_weight"] = composition.weight 
                    
                    for element in composition:
                        tmp_data_dict[element.value] = composition.get_atomic_fraction(element)

                    for material_element in self.material_elements:
                        if material_element not in list(tmp_data_dict.keys()):
                            tmp_data_dict[material_element] = 0 

                    self.data_dict[training_id] = tmp_data_dict

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
        #self.crystal_system_map = dict()
        #self.spacegroup_map = dict()

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
        #self.dataframe = self.dataframe[self.dataframe['band_gap'] != 0]

    def populate_data_dict_clean(self):
        """
        Populates data_dict_clean with training params only in self.non_element_keys.
        """
        for non_element_key in self.non_element_keys:
            self.data_dict_clean[non_element_key] = [ result[non_element_key] for (ID, result) in self.data_dict.items() ] 

        #if "crystal_system" in self.non_element_keys: 
        #    self.crystal_system_map = converters.create_non_numeric_map(self.data_dict_clean, "crystal_system")
        #    self.data_dict_clean["crystal_system"] = [self.crystal_system_map[value] for value in self.data_dict_clean["crystal_system"] ] 

        #if "spacegroup" in self.non_element_keys: 
        #    self.spacegroup_map = converters.create_non_numeric_map(self.data_dict_clean, "spacegroup")
        #    self.data_dict_clean["spacegroup"] = [self.spacegroup_map[value] for value in self.data_dict_clean["spacegroup"] ] 
    
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
        print("training params=",X_keys)

        X = np.asarray(self.dataframe[X_keys])
        y = np.asarray(self.dataframe['band_gap'])
                
        return train_test_split(X, y, test_size=test_size, shuffle= True) 
