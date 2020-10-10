import json
from os import listdir
from os.path import isfile, join

import sys
sys.path.insert(1, "../src/")
from converters import csv_to_json
from stoichiometry import get_norm_stoichiomertry

class BandGapDataset:
    def __init__(self,csv_path,json_path):
        self.csv_path = csv_path
        self.json_path = json_path
        self.data_dict = dict()
        self.data_IDs = list() 
        self.bad_formulas_IDs = list()
        self.convert_stored_csvs() 
        self.get_stored_data()
        self.set_stoichiometry()

    def convert_stored_csvs(self):
        csv_files = [f for f in listdir(self.csv_path) if isfile(join(self.csv_path, f))]
        json_files = [f for f in listdir(self.json_path) if isfile(join(self.json_path, f))]

        csv_files_stripped = [f.split(".csv")[0] for f in csv_files]
        json_files_stripped = [f.split(".json")[0] for f in json_files]

        for csv_file in csv_files_stripped:
            if csv_file not in json_files_stripped:
                csv_to_json(self.csv_path,self.json_path,fname=csv_file)
                print(f"created {csv_file}.json")

    def get_stored_data(self):
        json_files = [f for f in listdir(self.json_path) if isfile(join(self.json_path, f))]
        training_compounds = [f.split(".json")[0] for f in json_files]

        for training_compound in training_compounds:
            with open(f"{self.json_path}/{training_compound}.json") as fname:
                training_compound_results = dict(json.load(fname))
                for ID, training_compound_result in training_compound_results.items():
                    self.data_dict[ID] = training_compound_result

        self.data_IDs = list(self.data_dict.keys())

    def set_stoichiometry(self):
        for ID, result in self.data_dict.items():
            formula = result["formula"]
            if "(" not in formula:
                norm_stoichiometry = get_norm_stoichiomertry(formula)
                self.data_dict[ID]["stoichiometry"] = norm_stoichiometry
            else:
                self.bad_formulas_IDs.append(ID)

        for ID in self.bad_formulas_IDs:
            del self.data_dict[ID]
