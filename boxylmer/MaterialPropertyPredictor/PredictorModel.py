from abc import ABC, abstractmethod
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor


class AbstractBandGapModel(ABC):
    
    def __init__(self):
        pass
    
    @abstractmethod
    def train(self, material_data_loader):
        pass

    def predict(self, material_data_loader, test_data_only=True):
        if test_data_only:
            x = material_data_loader.get_test_data()
        else:
            x = material_data_loader.get_model_inputs()
        y = self._predict(x)
        return y

    @abstractmethod
    def _predict(self, x):
        pass

    def parity(self, material_data_loader, test_data_only=True):
        pred = self.predict(material_data_loader, test_data_only=test_data_only)
        if test_data_only:
            x, y = material_data_loader.get_test_data()
        else:
            x = material_data_loader.get_model_inputs()
            y = material_data_loader.get_model_outputs()

        return y, pred



class RandomForestBandGapModel(AbstractBandGapModel):
    def __init__(self):
        super.__init__()

    def train(self, material_data_loader):
        pass

    def _predict(self, x):
        pass


