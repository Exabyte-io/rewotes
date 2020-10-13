import pandas as pd
import numpy as np
from . import stoichiometry as stoichiometry 

import pymatgen as mg
from pymatgen.ext.matproj import MPRester

class Material:
    """
    Class used for user to create new material object for predictin of it's material.
    A formula is required input, and each additional input added will be used as a 
    training parameter.

    Args:
        formula (str)

    Kwargs:
        spacegroup (str)
        formation_energy (float)
        E_above_hull (float)
        volume (float)
        Nsites (int)
        density (float)
        crystal_system (str)
    """
    def __init__(self, formula, **kwargs): 
        self.formula = formula 
        self.features = dict(**kwargs)

class MaterialPredictionData:
    """
    Class used to create input data for each material. 
        - training parameters selected based on user input
        - material stoichiometry set 

    Arguments:
        material (Material object)
        symbols (list of str)
    """
    def __init__(self, material, symbols):
        self.symbols = symbols
        self.composition = mg.Composition(material.formula)
        self.molecular_weight = self.composition.weight
        self.material_stoichiometry = { element.value: self.composition.get_atomic_fraction(element) for element in self.composition }
        self.prediction_data = [ value for feature, value in material.features.items() ] 
        self.material_elements = list(self.material_stoichiometry.keys())   
        self.prediction_data.append(self.molecular_weight)
        for symbol in self.symbols:
            value = self.material_stoichiometry[symbol] if symbol in self.material_elements else 0
            self.prediction_data.append(value)
