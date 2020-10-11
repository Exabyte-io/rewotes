import numpy as np
from sklearn import linear_model, model_selection
import sys
sys.path.insert(1, "../src/")
from periodictable import PeriodicTable
from material import MaterialPredictionData
from trainingdata import BandGapDataset, BandGapDataFrame

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

class BandGapPrediction:
    """
    Arguments:
        csv_path (str)
        json_path (str)
        material (Material object)
        periodic_table (PeriodicTable object)
    """
    def __init__(self,csv_path,json_path,material,periodic_table):
        self.material = material
        self.symbols = periodic_table.symbols
        self.material_prediction_data = MaterialPredictionData(self.material, self.symbols)
        self.material_training_params = [ param for (param,value) in self.material.params.items() if value is not None ]
        self.band_gap_dataset = BandGapDataset(csv_path,json_path)  
        self.band_gap_dataframe_obj = BandGapDataFrame(self.band_gap_dataset.data_dict, self.symbols, self.material_training_params)
        self.make_prediction()

    def map_non_numeric_params(self):
        # mapping non-numeric params to numeric values which were created
        # during the creation of the training data
        non_numeric_params = {
            "crystal_system": self.band_gap_dataframe_obj.crystal_system_map,
            "spacegroup": self.band_gap_dataframe_obj.spacegroup_map,
        }
        for param, band_gap_map in non_numeric_params.items():
            if param in self.material_training_params:
                idx = self.material_training_params.index(param)
                value = self.material_prediction_data.prediction_data[idx]
                self.material_prediction_data.prediction_data[idx] = band_gap_map[value]            

    def train_model(self):
        X_train, X_test, y_train, y_test = self.band_gap_dataframe_obj.get_train_test_splits()

        # choose alpha
        # https://github.com/marcopeix/ISL-Ridge-Lasso/blob/master/Lasso%20and%20Ridge%20Regression.ipynb
        parameters = {
           'alpha': np.linspace(1e-3,1e1,50),
        }
        ridge_regressor = model_selection.GridSearchCV(linear_model.Ridge(), parameters, scoring='neg_mean_squared_error', cv=10)
        ridge_regressor.fit(X_train, y_train)
        alpha_choice = ridge_regressor.best_params_["alpha"]
        print(f"ridge regression alpha_choice = {alpha_choice}")

        model = linear_model.Ridge(alpha=alpha_choice)
        model.fit(X_train, y_train)

        self.X_train = X_train 
        self.X_test = X_test 
        self.y_train = y_train 
        self.y_test = y_test 
        self.model = model 
        self.model_score = model.score(X_test, y_test)
        self.model_weights = model.coef_

    def make_prediction(self):
        self.map_non_numeric_params()
        self.train_model()

        # predicting the bandgap
        this_prediction_data = np.asarray(self.material_prediction_data.prediction_data)
        this_prediction_data = np.reshape(this_prediction_data, (1,np.shape(this_prediction_data)[0]))

        self.predicted_band_gap = self.model.predict(this_prediction_data)[0]
