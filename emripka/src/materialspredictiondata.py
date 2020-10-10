import pandas as pd
import numpy as np

class MaterialsPredictionData:
    def __init__(self, materials_dataset, symbols):
        self.materials = materials_dataset.materials 
        self.materials_dict = materials_dataset.materials_dict 
        self.symbols = symbols
        self.prediction_data = dict()
        for material in self.materials:
            self.prediction_data[material.formula] = [ value for (param, value) in material.params.items() if value is not None ]
            self.populate_stoichiometry(material.formula)

    def populate_stoichiometry(self,formula):
        material_stoichiometry = self.materials_dict[formula]["stoichiometry"]  
        elements = list(material_stoichiometry.keys())   
        for symbol in self.symbols:
            value = material_stoichiometry[symbol] if symbol in elements else 0
            self.prediction_data[formula].append(value)
