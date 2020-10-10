import pandas as pd

class BandGapDataFrame:
    def __init__(self, data_dict, symbols):
        self.data_dict = data_dict
        self.symbols = symbols
        # create new dictionary
        self.data_dict_clean = {
            "ID": [ ID for (ID, result) in self.data_dict.items() ],
        }
        self.non_element_keys = ["band_gap__eV","density__gm_per_cc","volume","formation_energy__eV"] 
        self.populate_data_dict_clean()
        # create dataframe from data_dict_clean 
        self.dataframe = pd.DataFrame(self.data_dict_clean) 

    def populate_data_dict_clean(self):
        for non_element_key in self.non_element_keys:
            self.data_dict_clean[non_element_key] = [result[non_element_key] for (ID, result) in self.data_dict.items() ] 

        for symbol in self.symbols:
            self.data_dict_clean[symbol] = list() 

        self.populate_stoichiometry()
    
    def populate_stoichiometry(self):
        # populate the dictionary with keys of all symbols of elements
        for ID, result in self.data_dict.items():
            elements = list(result["stoichiometry"].keys())   
            for symbol in self.symbols:
                value = result["stoichiometry"][symbol] if symbol in elements else 0
                self.data_dict_clean[symbol].append(value)
