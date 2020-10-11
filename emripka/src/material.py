import pandas as pd
import numpy as np
import sys
sys.path.insert(1, "../src/")
from stoichiometry import get_norm_stoichiomertry

class Material:
    def __init__(self,formula,spacegroup=None,formation_energy=None,E_above_hull=None,
                    volume=None,Nsites=None,density=None,crystal_system=None):
        self.formula = formula
        self.spacegroup = spacegroup
        self.formation_energy = formation_energy
        self.E_above_hull = E_above_hull
        self.volume = volume
        self.Nsites = Nsites
        self.density = density
        self.crystal_system = crystal_system
        self.params = { 
            "spacegroup": self.spacegroup,
            "formation_energy__eV": self.formation_energy, 
            "E_above_hull__eV": self.E_above_hull,
            "volume": self.volume, 
            "Nsites": self.Nsites,
            "density__gm_per_cc": self.density,
            "crystal_system": self.crystal_system, 
        } 
        self.training_params = [ param for (param,value) in self.params.items() if value is not None ]

class MaterialPredictionData:
    """
    - the user creates an object of this type to use the package
    - takes in a MaterialsDatset type object which contains materials and their params to use for prediction
    - houses array of data to run through the model to predict the bandgap
    - also houses information needed to create the correct training data based on the user input
    """
    def __init__(self, material, symbols):
        self.symbols = symbols
        self.material_stoichiometry = get_norm_stoichiomertry(material.formula)  

        self.prediction_data = [ value for (param, value) in material.params.items() if value is not None ]

        self.elements = list(self.material_stoichiometry.keys())   
        for symbol in self.symbols:
            value = self.material_stoichiometry[symbol] if symbol in self.elements else 0
            self.prediction_data.append(value)
