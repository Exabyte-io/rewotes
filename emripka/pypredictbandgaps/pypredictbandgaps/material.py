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
        a (float)
        b (float)
        c (float)
        alpha (float)
        beta (float)
        gamma (float)
        volume (float)
    """
    def __init__(self, formula, **kwargs): 
        self.formula = formula 
        self.composition = mg.Composition(self.formula)
        self.features = dict(**kwargs)

class MaterialPredictionData:
    """
    Class used to create input data for each material. 
        - training parameters selected based on user input
        - material molecular weight and stoichiometry set 

    Args:
        material (Material)
    """
    def __init__(self, material):

        #self.molecular_weight = material.composition.weight
        self.material_stoichiometry = { element.value: material.composition.get_atomic_fraction(element) for element in material.composition }

        # stores numeric data for prediction
        self.prediction_data = [ value for feature, value in material.features.items() ] 

        # list of parameters which are trained on/used for prediciton
        self.training_params = list(material.features.keys())

        self.material_elements = list(self.material_stoichiometry.keys())   
        for element in self.material_elements:
            self.training_params.append(element)

        for element in material.composition:
            self.prediction_data.append(material.composition.get_atomic_fraction(element))
        #self.prediction_data.append(self.molecular_weight)
