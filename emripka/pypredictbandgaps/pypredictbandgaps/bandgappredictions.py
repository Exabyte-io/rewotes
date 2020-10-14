import numpy as np
from sklearn import linear_model, model_selection
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

from . import periodictable as periodictable
from . import material as material_module
from . import trainingdata as trainingdata

import json
import shutil
from datetime import datetime
import os
this_dir, this_filename = os.path.split(__file__)


class BandGapPredictions:
    """ 
    Class used to predict the band gap of a set of materials. A unique model is 
    trained for each material, as each material can have unique input parameters.
    
    Args:
        materials (list of Material objects)
        training_data (list of str)
    """
    def __init__(self, materials, model_type, training_data=None):
        self.band_gap_prediction_objects = { material.formula: BandGapPrediction(material, model_type) for material in materials }

    def move_user_data(self):
        """ 
        Moves the user_data.json file to /data/training/user_data_json/ and rename to the date and time
        of the prediciton after predictions
        """
        json_path = this_dir+"/data/training/materialsproject_json"
        user_data_path = this_dir+"/data/training/user_data_json"
        fname = create_user_data_filename()   
        shutil.move(f"{json_path}/user_data.json", f"{user_data_path}/{fname}.json")

        # create empty user_data.json file
        with open(f"{json_path}/user_data.json","w") as fname:
            json.dump(dict(),fname,indent=4)
        fname.close()

class BandGapPrediction:
    """
    Unique training dataset created and used to train the model for a user-defined 
    Material object. 

    Args:
        material (Material)
    """
    def __init__(self, material, model_type):
        self.material = material
        self.model_type = model_type
        self.periodic_table = periodictable.PeriodicTable()

        self.material_prediction_data = material_module.MaterialPredictionData(self.material)
        self.material_training_params = self.material_prediction_data.training_params

        self.band_gap_dataset_obj = trainingdata.BandGapDataset(self.material)  
        self.band_gap_dataframe_obj = trainingdata.BandGapDataFrame(self.band_gap_dataset_obj.data_dict, self.periodic_table.symbols, self.material_training_params)

        #self.map_non_numeric_params()
        self.train_model()

    def map_non_numeric_params(self):
        """
        Maps non-numeric params to numeric values which were created
        during the creation of the training data.
        """
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
        if self.model_type == "ridge_regression":
            self.train_model_ridge_regression()
        elif self.model_type == "svm":
            self.train_model_svm()

    def choose_alpha(self, X_train, y_train):
        """
        Alpha parameter selection for Ridge Regression.
        https://github.com/marcopeix/ISL-Ridge-Lasso/blob/master/Lasso%20and%20Ridge%20Regression.ipynb

        Args:
            X_train (arr)
            y_train (arr)
        """
        cv = 10 if np.shape(X_train)[0] > 10 else np.shape(X_train)[0] 
        parameters = {
           'alpha': np.linspace(1e-3,1e1,50),
        }
        ridge_regressor = model_selection.GridSearchCV(linear_model.Ridge(), parameters, scoring='neg_mean_squared_error', cv=cv)
        ridge_regressor.fit(X_train, y_train)
        alpha_choice = ridge_regressor.best_params_["alpha"]

        print(f"ridge regression alpha_choice = {alpha_choice}")
        return alpha_choice

    def train_model_ridge_regression(self):
        """
        Ridge Regression model training.   
        """
        print("Training the model...")
        X_train, X_test, y_train, y_test = self.band_gap_dataframe_obj.get_train_test_splits()

        alpha_choice = self.choose_alpha(X_train, y_train)
        model = linear_model.Ridge(alpha=alpha_choice)
        model.fit(X_train, y_train)

        self.X_train = X_train 
        self.X_test = X_test 
        self.y_train = y_train 
        self.y_test = y_test 
        self.model = model 
        self.model_score = model.score(X_test, y_test)
        self.model_weights = model.coef_

        print("Making predictions...")
        this_prediction_data = np.asarray(self.material_prediction_data.prediction_data)
        this_prediction_data = np.reshape(this_prediction_data, (1,np.shape(this_prediction_data)[0]))
        print("prediciton data:",this_prediction_data)
        self.predicted_band_gap = self.model.predict(this_prediction_data)[0]
        print("predicted bandgap:", self.predicted_band_gap) 

    def train_model_svm(self):
        """
        Support Vector Machine model training.   
        """
        X, y = self.band_gap_dataframe_obj.get_all_train_data()
        sc_X = StandardScaler()
        sc_y = StandardScaler()

        y = np.reshape(y, (np.shape(y)[0],1))

        X = sc_X.fit_transform(X)
        y = sc_y.fit_transform(y)

        self.X = X
        self.y = y 

        regressor = SVR(kernel='rbf')
        model = regressor.fit(X, y)
        self.model = model

        print("Making predictions...")
        this_prediction_data = np.asarray(self.material_prediction_data.prediction_data)
        this_prediction_data = np.reshape(this_prediction_data, (1,np.shape(this_prediction_data)[0]))
        print("prediciton data:",this_prediction_data)
        self.predicted_band_gap = sc_y.inverse_transform((regressor.predict (sc_X.transform(this_prediction_data))))
        print("predicted bandgap:", self.predicted_band_gap) 
