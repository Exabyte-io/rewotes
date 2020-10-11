import pandas as pd
import sys
sys.path.insert(1, "../src/")
from converters import create_non_numeric_map

class BandGapDataFrame:
    def __init__(self, data_dict, symbols,material_training_params):
        self.data_dict = data_dict
        self.symbols = symbols
        self.crystal_system_map = dict()
        self.spacegroup_map = dict()
        # create new dictionary
        self.data_dict_clean = {
            "ID": [ ID for (ID, result) in self.data_dict.items() ],
        }
        self.non_element_keys = ["band_gap__eV",]
        for param in material_training_params:
            self.non_element_keys.append(param)
        self.populate_data_dict_clean()

        # create dataframe from data_dict_clean 
        self.dataframe = pd.DataFrame(self.data_dict_clean) 
        # dropping rows which have a bandgap of zero
        #self.dataframe = self.dataframe[self.dataframe['band_gap__eV'] != 0]

    def populate_data_dict_clean(self):
        for non_element_key in self.non_element_keys:
            self.data_dict_clean[non_element_key] = [result[non_element_key] for (ID, result) in self.data_dict.items() ] 

        for symbol in self.symbols:
            self.data_dict_clean[symbol] = list() 

        self.populate_stoichiometry()

        if "crystal_system" in self.non_element_keys: 
            self.crystal_system_map = create_non_numeric_map(self.data_dict_clean, "crystal_system")
            self.data_dict_clean["crystal_system"] = [self.crystal_system_map[value] for value in self.data_dict_clean["crystal_system"] ] 

        if "spacegroup" in self.non_element_keys: 
            self.spacegroup_map = create_non_numeric_map(self.data_dict_clean, "spacegroup")
            self.data_dict_clean["spacegroup"] = [self.spacegroup_map[value] for value in self.data_dict_clean["spacegroup"] ] 
    
    def populate_stoichiometry(self):
        # populate the dictionary with keys of all symbols of elements
        for ID, result in self.data_dict.items():
            elements = list(result["stoichiometry"].keys())   
            for symbol in self.symbols:
                value = result["stoichiometry"][symbol] if symbol in elements else 0
                self.data_dict_clean[symbol].append(value)
