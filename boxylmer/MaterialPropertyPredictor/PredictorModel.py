from abc import ABC, abstractmethod
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import RandomizedSearchCV


class AbstractBandGapModel(ABC):
    
    def __init__(self):
        """
        Abstract band gap prediction model.
        
        Defines the interface and some common functionality for Band Gap prediction.

        Accepts a generic MaterialDataLoader, which will predict using however many features that data loader provides.
        """
        self.model = None
        pass
    
    @abstractmethod
    def fit(self, material_data_loader):
        """
        Fit the model to the training data in the data loader.
        Optionally, a printed progress bar may be supplied.

        Args:
            material_data_loader (AbstractDataLoader): Data loader object containing training and testing data.
        """
        pass

    def predict(self, material_data_loader, test_data_only=True):
        """
        Predict band gap based on existing test data within the loader. Optionally predict using ALL data if test_data_only is false. 
        - By default, a seed is used, so data should be the same at every call. 

        Args:
            material_data_loader (AbstractDataLoader): Data loader object containing material data.
            test_data_only (bool): If True (default), use only test data for predictions.

        Returns:
            numpy.ndarray: Predicted band gap values.
        """
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
        """
        Get predicted and actual band gap values for the data loader's test data or all data.
        
        Args:
            material_data_loader (AbstractDataLoader): data loader object.
            test_data_only (bool): If True, use only previously unseen test data.
            
        Returns:
            tuple: (Actual band gap values, predicted band gap values)
        """
        if test_data_only:
            x, y = material_data_loader.get_test_data()
        else:
            x = material_data_loader.get_model_inputs()
            y = material_data_loader.get_model_outputs()

        pred = self._predict(x)

        return y, pred
    
    @abstractmethod
    def fit_hyperparameters(self, material_data_loader):
        """
        Fit hyperparameters of the model using training data, then train the model on the optimized hyperparameters, which are also printed to console. 
        
        Args:
            material_data_loader (AbstractDataLoader): Data loader object.

        Returns: 
            None
        """
        pass

class RandomForestBandGapModel(AbstractBandGapModel):
    """Predict band gap using random forest regression."""
    def __init__(self, **kwargs):
        super().__init__()
        self.model = RandomForestRegressor(
            n_estimators=160, 
            min_samples_split=4, 
            min_samples_leaf=1, 
            max_depth=None, 
            bootstrap=True,
            **kwargs
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
    """Predict band gap using gradient boosted trees."""
    def __init__(self, **kwargs):
        super().__init__()
        self.model = GradientBoostingRegressor(
            n_estimators=180, 
            min_samples_split=3, 
            min_samples_leaf=8, 
            max_depth=3, 
            learning_rate=0.1,
            **kwargs,
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
        
        return None