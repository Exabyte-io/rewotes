from abc import ABC, abstractmethod
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import RandomizedSearchCV


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
    
    @abstractmethod
    def fit_hyperparameters(self, material_data_loader):
        pass

class RandomForestBandGapModel(AbstractBandGapModel):
    def __init__(self):
        super().__init__()
        self.model = RandomForestRegressor(
            n_estimators=50, 
            min_samples_split=4, 
            min_samples_leaf=2, 
            max_depth=30, 
            bootstrap=True
            )

        
    def fit(self, material_data_loader):
        train_x, train_y = material_data_loader.get_train_data()
        self.model.fit(train_x, train_y)
        return None

    def _predict(self, x):
        return self.model.predict(x)
    

    def fit_hyperparameters(self, material_data_loader, n_iter=50):
        param_dist = {
            'n_estimators': np.arange(10, 200, 10),
            'max_depth': [None] + list(np.arange(5, 50, 5)),
            'min_samples_split': np.arange(2, 10, 1),
            'min_samples_leaf': np.arange(1, 10, 1),
            'bootstrap': [True, False]
        }
        
        random_search = RandomizedSearchCV(
            estimator=self.model,
            param_distributions=param_dist,
            n_iter=n_iter,
            scoring='neg_mean_squared_error', # really annoying convention, I'm going to forget this later and be upset when I do it again and it moves **away** from local minima
            verbose=1,
            n_jobs=-1
        )
        
        train_x, train_y = material_data_loader.get_train_data()
        
        random_search.fit(train_x, train_y) # very confused initially but the docs say it does a validation split internally, so we're OK

        self.model = random_search.best_estimator_

        print(f"Best hyperparameters found: {random_search.best_params_} at RMSE = {-random_search.best_score_}")

# if this performs well, we would want to replace it with the more performant XGBoost model, 
#   which is fundamentally the same but can be parallelized easily and may even have better optimization methods
# more info: https://stats.stackexchange.com/questions/282459/xgboost-vs-python-sklearn-gradient-boosted-trees
class GradientBoostingBandGapModel(AbstractBandGapModel):
    def __init__(self):
        super().__init__()
        self.model = GradientBoostingRegressor(
            n_estimators=170, 
            min_samples_split=4, 
            min_samples_leaf=6, 
            max_depth=4, 
            learning_rate=0.1
        )
        

    def fit(self, material_data_loader):
        train_x, train_y = material_data_loader.get_train_data()
        self.model.fit(train_x, train_y)
        return None

    def _predict(self, x):
        return self.model.predict(x)
    
    def fit_hyperparameters(self, material_data_loader, n_iter=50):
        param_dist = {
            'n_estimators': np.arange(10, 200, 10),
            'max_depth': [None] + list(np.arange(3, 16, 1)),
            'min_samples_split': np.arange(2, 10, 1),
            'min_samples_leaf': np.arange(1, 10, 1),
            'learning_rate': [0.01, 0.05, 0.1, 0.2, 0.5]
        }
        
        random_search = RandomizedSearchCV(
            estimator=self.model,
            param_distributions=param_dist,
            n_iter=n_iter,
            scoring='neg_mean_squared_error', 
            verbose=1,
            n_jobs=-1
        )
        
        train_x, train_y = material_data_loader.get_train_data()
        
        random_search.fit(train_x, train_y) 

        self.model = random_search.best_estimator_

        print(f"Best hyperparameters found: {random_search.best_params_} at RMSE = {-random_search.best_score_}")