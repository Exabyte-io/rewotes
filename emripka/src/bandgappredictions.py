import numpy as np
from sklearn import linear_model, model_selection, metrics
import sys
sys.path.insert(1, "../src/")
from bandgapdataset import BandGapDataset
from bandgapdataframe import BandGapDataFrame

class BandGapPredictions:
    # To-do: create a class for each material to use, which are housed inside of BandGapPredictions
    def __init__(self,materials_to_predict,csv_path,json_path,symbols,materials_dataset,materials_prediction_data):
        self.materials_to_predict = materials_to_predict
        self.symbols = symbols
        self.predictions = { material: dict() for material in self.materials_to_predict }
        self.materials_dataset = materials_dataset
        self.materials_prediction_data = materials_prediction_data
        self.materials_training_params = dict()
        self.band_gap_dataset = BandGapDataset(csv_path,json_path)  
        for material in materials_to_predict:
            self.make_prediction(material)

    def make_band_gap_dataframe(self, material):
        self.materials_training_params[material] = [ 
            param for (param,value) in self.materials_dataset.materials_dict[material]["params"].items() if value is not None
        ]

        band_gap_dataframe_obj = BandGapDataFrame(self.band_gap_dataset.data_dict, self.symbols, self.materials_training_params[material])

        non_element_keys = band_gap_dataframe_obj.non_element_keys 
        self.predictions[material]["params"] = non_element_keys[1:]
        
        return band_gap_dataframe_obj

    def map_non_numeric_params(self, material, band_gap_dataframe_obj):
        # mapping non-numeric params to numeric values which were created
        # during the creation of the training data
        non_numeric_params = {
            "crystal_system": band_gap_dataframe_obj.crystal_system_map,
            "spacegroup":band_gap_dataframe_obj.spacegroup_map,
        }
        for param, band_gap_map in non_numeric_params.items():
            if param in self.materials_training_params[material]:
                idx = self.materials_training_params[material].index(param)
                value = self.materials_prediction_data.prediction_data[material][idx]
                self.materials_prediction_data.prediction_data[material][idx] = band_gap_map[value]            

    def train_model(self, material, band_gap_dataframe_obj):
        X_train, X_test, y_train, y_test = band_gap_dataframe_obj.get_train_test_splits()
        model = linear_model.Ridge(alpha=0.5)
        model.fit(X_train, y_train)

        self.predictions[material]["X_train"] = X_train 
        self.predictions[material]["X_test"] = X_test 
        self.predictions[material]["y_train"] = y_train 
        self.predictions[material]["y_test"] = y_test 
        self.predictions[material]["model"] = model 
        self.predictions[material]["model_score"] = model.score(X_test, y_test)
        self.predictions[material]["model_weights"] = model.coef_

        return model

    def make_prediction(self, material):
        band_gap_dataframe_obj = self.make_band_gap_dataframe(material)
        self.map_non_numeric_params(material, band_gap_dataframe_obj)
        trained_model = self.train_model(material, band_gap_dataframe_obj)

        # predicting the bandgap
        this_prediction_data = np.asarray(self.materials_prediction_data.prediction_data[material])
        this_prediction_data = np.reshape(this_prediction_data, (1,np.shape(this_prediction_data)[0]))

        self.predictions[material]["band_gap__eV"] = trained_model.predict(this_prediction_data)[0]
