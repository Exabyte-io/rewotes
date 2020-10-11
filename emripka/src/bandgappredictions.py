import numpy as np
from sklearn import linear_model, model_selection, metrics
import sys
sys.path.insert(1, "../src/")
from periodictable import PeriodicTable
from bandgapprediction import BandGapPrediction

class BandGapPredictions:
    """
    - prediction made for each material
    - unique model is trained for each material, as each material can have unique input parameters
    """
    def __init__(self,csv_path,json_path,materials_list):
        self.periodic_table = PeriodicTable()
        self.band_gap_prediction_objects = { 
            material.formula: BandGapPrediction(csv_path,json_path,material,self.periodic_table) for material in materials_list 
        }
