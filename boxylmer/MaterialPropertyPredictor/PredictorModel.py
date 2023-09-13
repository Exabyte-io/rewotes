from abc import ABC, abstractmethod
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor


class AbstractBandGapModel(ABC):
    
    def __init__(self):
        self.model = None
        pass
    
    @abstractmethod
    def fit(self, material_data_loader):
        pass

    def predict(self, material_data_loader, test_data_only=True):
        if test_data_only:
            x, _ = material_data_loader.get_test_data()
        else:
            x = material_data_loader.get_model_inputs()
        y = self._predict(x)
        return y

    @abstractmethod
    def _predict(self, x):
        pass

    def parity(self, material_data_loader, test_data_only=True):
        if test_data_only:
            x, y = material_data_loader.get_test_data()
        else:
            x = material_data_loader.get_model_inputs()
            y = material_data_loader.get_model_outputs()

        pred = self._predict(x)

        return y, pred

class RandomForestBandGapModel(AbstractBandGapModel):
    def __init__(self):
        super().__init__()
        self.model = RandomForestRegressor(n_estimators=150)
        

    def fit(self, material_data_loader):
        train_x, train_y = material_data_loader.get_train_data()
        self.model.fit(train_x, train_y)
        return None

    def _predict(self, x):
        return self.model.predict(x)



class GradientBoostingBandGapModel(AbstractBandGapModel):
    def __init__(self):
        super().__init__()
        self.model = GradientBoostingRegressor(n_estimators=150, learning_rate=0.1)
        

    def fit(self, material_data_loader):
        train_x, train_y = material_data_loader.get_train_data()
        self.model.fit(train_x, train_y)
        return None

    def _predict(self, x):
        return self.model.predict(x)
    