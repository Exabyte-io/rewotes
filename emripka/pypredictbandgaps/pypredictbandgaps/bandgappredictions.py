import numpy as np
from sklearn import linear_model, model_selection
from . import periodictable as periodictable
from . import material as material_module
from . import trainingdata as trainingdata

class BandGapPredictions:
    """ 
    Class used to predict the band gap of a set of materials. A unique model is 
    trained for each material, as each material can have unique input parameters.
    
    Arguments:
        materials_list (list of Material objects):
            material (Material object)
        input_training_data (list of TrainingData objects)
        use_database_data (bool): decides useage of the package
            True: the database data is used to train the model, along with any input TrainingData objects 
            False: the database data is not used to train the model, and only input TrainingData is
    To-do:
    """
    def __init__(self, materials_list,input_training_data=None,use_database_data=True):
        self.periodic_table = periodictable.PeriodicTable()
        self.use_database_data = use_database_data
        if input_training_data is not None:
            for input_data in input_training_data:
                input_data.store_data()
        else:
            self.use_database_data = True
        self.band_gap_prediction_objects = { 
            material.formula: BandGapPrediction(material, self.periodic_table, use_database_data) for material in materials_list 
        }

class BandGapPrediction:
    """
    Unique training dataset created and used to train the model for a user-defined 
    Material object. 

    Arguments:
        material (Material object)
        periodic_table (PeriodicTable object)
        use_database_data (bool): decides useage of the package
            True: the database data is used to train the model, along with any input TrainingData objects 
            False: the database data is not used to train the model, and only input TrainingData is
    """
    def __init__(self, material, periodic_table, use_database_data):
        self.material = material
        self.material_prediction_data = material_module.MaterialPredictionData(self.material, periodic_table.symbols, periodic_table)
        self.material_training_params = [ param for (param,value) in self.material.params.items() if value is not None ]
        self.band_gap_dataset_obj = trainingdata.BandGapDataset(periodic_table, use_database_data)  
        self.band_gap_dataframe_obj = trainingdata.BandGapDataFrame(self.band_gap_dataset_obj.data_dict, periodic_table.symbols, self.material_training_params)
        self.map_non_numeric_params()
        self.train_model()
        self.make_prediction()

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
        """
        Alpha parameter selection and Ridge Regression model training.   
        """
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
        """
        Reshapes the user-input material data and makes a prediction of the 
        band gap of the material using the trained model.
        """
        # predicting the bandgap
        this_prediction_data = np.asarray(self.material_prediction_data.prediction_data)
        this_prediction_data = np.reshape(this_prediction_data, (1,np.shape(this_prediction_data)[0]))

        self.predicted_band_gap = self.model.predict(this_prediction_data)[0]
