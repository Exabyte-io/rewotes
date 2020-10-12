import pandas as pd
import numpy as np
from . import stoichiometry as stoichiometry 

class Material:
    """
    Class used for user to create new material object for predictin of it's material.
    A formula is required input, and each additional input added will be used as a 
    training parameter.

    Arguments:
        formula (str)
        spacegroup (str)
        formation_energy (float)
        E_above_hull (float)
        volume (float)
        Nsites (int)
        density (float)
        crystal_system (str)
    """
    def __init__(self,formula,spacegroup=None,formation_energy=None,E_above_hull=None,
                    volume=None,Nsites=None,density=None,crystal_system=None):
        self.formula = formula
        self.params = { 
            "spacegroup": spacegroup,
            "formation_energy__eV": formation_energy, 
            "E_above_hull__eV": E_above_hull,
            "volume": volume, 
            "Nsites": Nsites,
            "density__gm_per_cc": density,
            "crystal_system": crystal_system, 
        } 
        self.training_params = [ param for (param,value) in self.params.items() if value is not None ]

class MaterialPredictionData:
    """
    Class used to create input data for each material. 
        - training parameters selected based on user input
        - material stoichiometry set 

    Arguments:
        material (Material object)
        symbols (list of str)
        periodic_table (PeriodicTable object)
    """
    def __init__(self, material, symbols, periodic_table):
        self.symbols = symbols
        self.molecular_weight = stoichiometry.get_molecular_weight(material.formula, periodic_table)
        self.material_stoichiometry = stoichiometry.get_norm_stoichiomertry(material.formula)  
        self.prediction_data = [ value for (param, value) in material.params.items() if value is not None ]
        self.material_elements = list(self.material_stoichiometry.keys())   
        self.prediction_data.append(self.molecular_weight)
        for symbol in self.symbols:
            value = self.material_stoichiometry[symbol] if symbol in self.material_elements else 0
            self.prediction_data.append(value)
