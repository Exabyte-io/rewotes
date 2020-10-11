import sys
sys.path.insert(1, "../src/")
from stoichiometry import get_norm_stoichiomertry
from converters import create_non_numeric_map

class MaterialsDataset:
    def __init__(self,materials):
        self.materials = materials
        self.materials_dict = { material.formula: dict() for material in self.materials }
        for material in materials:
            formula = material.formula
            self.materials_dict[formula]["params"] = material.params 
            self.materials_dict[formula]["stoichiometry"] = get_norm_stoichiomertry(formula)
