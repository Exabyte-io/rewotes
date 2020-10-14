import pandas as pd
import numpy as np
from . import stoichiometry as stoichiometry 

import pymatgen as mg
from pymatgen.symmetry.groups import sg_symbol_from_int_number

class Material:
    """
    Class used for user to create new material object for predictin of it's material.
    A formula is required input, and each additional input added will be used as a 
    training parameter.

    Args:
        formula (str)

    Kwargs:
        spacegroup (str)
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
        self.spacegroup = self.features["spacegroup"]
        spacegroup_map = { sg_symbol_from_int_number(ii):ii for ii in range(1,231)}
        self.features["spacegroup"] = spacegroup_map[self.features["spacegroup"]]

class MaterialPredictionData:
    """
    Class used to create input data for each material. 
        - training parameters selected based on user input
        - material molecular weight and stoichiometry set 

    Args:
        material (Material)
        model_type (str): choose from the following models:
            ["ridge_regression", "svm", "decision_tree", "random_forest"]
    """
    def __init__(self, material, model_type):
        self.molecular_weight = material.composition.weight
        self.material = material
        self.material_stoichiometry = { element.value: material.composition.get_atomic_fraction(element) for element in material.composition }
        self.material_features = material.features
        self.model_type = model_type
        self.assign_params()

    def assign_params(self):
        """
        Chooses model params based on user choice.
        """
        if self.model_type == "ridge_regression":
            self.assign_params_ridge_regression()
        elif self.model_type == "svm":
            self.assign_params_svm()
        elif self.model_type == "decision_tree":
            self.assign_params_decision_tree()
        elif self.model_type == "random_forest":
            self.assign_params_random_forest()

    def assign_params_ridge_regression(self):
        """
        Parameter selection for Ridge Regression.
        """
        #for param in ["a","b","c","alpha","beta","gamma"]:
        for param in ["alpha","beta","gamma"]:
            del self.material_features[param]

        del self.material_features["spacegroup"]

        # stores numeric data for prediction
        self.prediction_data = [ value for feature, value in self.material_features.items() ] 

        # list of parameters which are trained on/used for prediciton
        self.training_params = list(self.material_features.keys())

        self.material_elements = list(self.material_stoichiometry.keys())   
        for element in self.material_elements:
            self.training_params.append(element)          # i.e. "Si"
            #self.training_params.append(element+"_group") # i.e "Si_group"
            #self.training_params.append(element+"_electronegativity") # i.e "Si_electronegativity"

        for element in self.material.composition:
            self.prediction_data.append(self.material.composition.get_atomic_fraction(element))
            tmp_element = mg.Element(element.value)
            #self.prediction_data.append(tmp_element.group)
            #self.prediction_data.append(tmp_element.X)

        #self.training_params.append("molecular_weight")
        #self.prediction_data.append(molecular_weight)

    def assign_params_random_forest(self):
        """
        Parameter selection for Random Forest.
        """
        #for param in ["a","b","c","alpha","beta","gamma"]:
        #    del self.material_features[param]

        # stores numeric data for prediction
        self.prediction_data = [ value for feature, value in self.material_features.items() ] 

        # list of parameters which are trained on/used for prediciton
        self.training_params = list(self.material_features.keys())

        self.material_elements = list(self.material_stoichiometry.keys())   
        for element in self.material_elements:
            self.training_params.append(element)          # i.e. "Si"
            self.training_params.append(element+"_group") # i.e "Si_group"
            self.training_params.append(element+"_electronegativity") # i.e "Si_electronegativity"

        for element in self.material.composition:
            self.prediction_data.append(self.material.composition.get_atomic_fraction(element))
            tmp_element = mg.Element(element.value)
            self.prediction_data.append(tmp_element.group)
            self.prediction_data.append(tmp_element.X)

        self.training_params.append("molecular_weight")
        self.prediction_data.append(self.molecular_weight)

    def assign_params_svm(self):
        """
        Parameter selection for SVM.
        """
        #for param in ["a","b","c","alpha","beta","gamma"]:
        #    del self.material_features[param]

        # stores numeric data for prediction
        self.prediction_data = [ value for feature, value in self.material_features.items() ] 

        # list of parameters which are trained on/used for prediciton
        self.training_params = list(self.material_features.keys())

        self.material_elements = list(self.material_stoichiometry.keys())   
        for element in self.material_elements:
            self.training_params.append(element)          # i.e. "Si"
            self.training_params.append(element+"_group") # i.e "Si_group"
            self.training_params.append(element+"_electronegativity") # i.e "Si_electronegativity"

        for element in self.material.composition:
            self.prediction_data.append(self.material.composition.get_atomic_fraction(element))
            tmp_element = mg.Element(element.value)
            self.prediction_data.append(tmp_element.group)
            self.prediction_data.append(tmp_element.X)

        self.training_params.append("molecular_weight")
        self.prediction_data.append(self.molecular_weight)

    def assign_params_decision_tree(self):
        """
        Parameter selection for Decision Tree.
        """
        #for param in ["a","b","c","alpha","beta","gamma"]:
        #    del self.material_features[param]

        # stores numeric data for prediction
        self.prediction_data = [ value for feature, value in self.material_features.items() ] 

        # list of parameters which are trained on/used for prediciton
        self.training_params = list(self.material_features.keys())

        self.material_elements = list(self.material_stoichiometry.keys())   
        for element in self.material_elements:
            self.training_params.append(element)          # i.e. "Si"
            self.training_params.append(element+"_group") # i.e "Si_group"
            self.training_params.append(element+"_electronegativity") # i.e "Si_electronegativity"

        for element in self.material.composition:
            self.prediction_data.append(self.material.composition.get_atomic_fraction(element))
            tmp_element = mg.Element(element.value)
            self.prediction_data.append(tmp_element.group)
            self.prediction_data.append(tmp_element.X)

        self.training_params.append("molecular_weight")
        self.prediction_data.append(self.molecular_weight)
