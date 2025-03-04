"""
This module contains classes for training and evaluating materials models.

The module includes the following classes:
- MaterialsModels: A base class for training and evaluating machine learning
models on materials data.
- XGBoostModels: A class for training and evaluating XGBoost models on
materials data.
- NeuralNetModels: A class for training and evaluating PyTorch models on
materials data.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import uniform, randint
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold, RandomizedSearchCV
from skorch import NeuralNetRegressor
import torch
from torch import nn
import xgboost as xgb


class MaterialsModels:
    """A class for training and evaluating materials models.

    This class provides a template for training and evaluating machine learning
    models on materials data. It includes methods for training models with
    cross-validation and randomized search for hyperparameter tuning, as well
    as evaluating the trained models.

    Args:
        x_train (np.ndarray): The training features.
        y_train (np.ndarray): The training labels.
        x_test (np.ndarray): The testing features.
        y_test (np.ndarray): The testing labels.
        save (bool, optional): Whether to save the trained model and metrics.
        Defaults to True.

    Attributes:
        x_train (np.ndarray): The training features.
        y_train (np.ndarray): The training labels.
        x_test (np.ndarray): The testing features.
        y_test (np.ndarray): The testing labels.
        save (bool): Whether to save the trained model and metrics.
        metrics (dict): The evaluation metrics.
        model: The trained model.
        _dir_out (Path): The output directory for saving the model and metrics.

    Raises:
        NotImplementedError: If the train_models method is not implemented.

    """

    def __init__(
        self,
        x_train: np.ndarray,
        y_train: np.ndarray,
        x_test: np.ndarray,
        y_test: np.ndarray,
        save: bool = True,
    ) -> None:
        """Initialize the MaterialsModels class.

        Args:
            x_train (np.ndarray): The training features.
            y_train (np.ndarray): The training labels.
            x_test (np.ndarray): The testing features.
            y_test (np.ndarray): The testing labels.
            save (bool, optional): Whether to save the trained model and
            metrics. Defaults to True.
        """

        self.x_train = x_train
        self.y_train = y_train
        self.x_test = x_test
        self.y_test = y_test
        self.save = save

        self.metrics = None
        self.model = None

        self._dir_out = Path("./models")

    def train_models(
        self, param_grid: dict = None, cv: int = 10, seed: int = 42
    ) -> dict:
        """Train XGBoost models with cross-validation and grid search.

        Args:
            param_grid (dict): The parameter grid for the grid search.
            Defaults to None.
            cv (int, optional): The number of cross-validation folds.
            Defaults to 5.
            seed (int, optional): The random seed for the train-test split.
            Defaults to 42.

        Returns:
            dict: The trained models.

        Raises:
            NotImplementedError: If the method is not implemented.
        """
        raise NotImplementedError("Method not implemented")

    def evaluate_model(self) -> dict:
        """Evaluate the trained model.

        Returns:
            dict: The evaluation metrics.
        """
        y_pred = self.model.predict(self.x_test)
        mse = mean_squared_error(self.y_test, y_pred)

        self.metrics = {
            "mean_squared_error": mse,
            "root_mean_squared_error": np.sqrt(mse),
        }
        return self.metrics


class XGBoostModels(MaterialsModels):
    """
    This class represents a set of XGBoost models for materials data.

    It inherits from the MaterialsModels base class and implements the
    train_models method specifically for XGBoost models. It uses
    cross-validation and grid search for training.

    Attributes:
        model: The trained model.
        metrics: The evaluation metrics for the model.
        x_train: The training data.
        y_train: The training labels.
        x_test: The test data.
        y_test: The test labels.
    """

    def train_models(
        self, param_grid: dict = None, cv: int = 5, seed: int = 42
    ) -> xgb.XGBRegressor:
        """Train XGBoost models with cross-validation and grid search.

        Args:
            param_grid (dict): The parameter grid for the grid search.
            Defaults to None.
            cv (int, optional): The number of cross-validation folds.
            Defaults to 5.
            seed (int, optional): The random seed for the train-test split.
            Defaults to 42.

        Returns:
            xgb.XGBRegressor: The best trained XGBoost model.
        """
        if param_grid is None:
            param_grid = {
                "n_estimators": randint(100, 1000),
                "max_depth": randint(3, 10),
                "subsample": uniform(0.5, 0.5),
                "colsample_bytree": uniform(0.5, 0.5),
            }

        model = xgb.XGBRegressor(
            n_estimators=1000,
            max_depth=7,
            learning_rate=0.1,
            colsample_bytree=0.8,
            subsample=0.8,
            n_jobs=1,
        )
        kfold = KFold(n_splits=cv, random_state=seed, shuffle=True)
        random_search = RandomizedSearchCV(
            model,
            param_distributions=param_grid,
            n_iter=20,
            scoring="neg_mean_squared_error",
            n_jobs=-1,
            cv=kfold.split(self.x_train, self.y_train),
            refit=True,
            verbose=3,
            random_state=seed,
        )
        random_search.fit(self.x_train, self.y_train)
        best_model = random_search.best_estimator_
        self.model = best_model

        params = random_search.best_params_
        print(f"Best Model:\n{best_model}")

        y_pred = best_model.predict(self.x_test)
        mse = mean_squared_error(self.y_test, y_pred)
        print(f"Mean Squared Error: {mse}")

        if self.save:
            self._dir_out.mkdir(exist_ok=True)
            best_model.save_model(self._dir_out / "xgboost_model.json")
            pd.DataFrame(params, index=[0]).to_csv(
                self._dir_out / "xgboost_params.csv", index=False
            )

        return best_model

    def evaluate_model(self) -> dict:
        """Evaluate the trained model.

        Returns:
            dict: The evaluation metrics.
        """
        super().evaluate_model()

        if self.save:
            self._dir_out.mkdir(exist_ok=True)
            pd.DataFrame(self.metrics, index=[0]).to_csv(
                self._dir_out / "xgboost_metrics.csv", index=False
            )

        return self.metrics


class NeuralNetModels(MaterialsModels):
    """
    This class represents a set of PyTorch models for materials data.

    It inherits from the MaterialsModels base class and implements the
    train_models method specifically for PyTorch models. It uses
    cross-validation and grid search for training.

    Attributes:
        model: The trained model.
        metrics: The evaluation metrics for the model.
        x_train: The training data.
        y_train: The training labels.
        x_test: The test data.
        y_test: The test labels.
    """

    class Net(nn.Module):
        """A simple feedforward neural network.

        Args:
            nn (Module): The PyTorch neural network module.
        """

        def __init__(
            self,
            input_size,
            hidden_size,
            output_size,
            num_layers,
            dropout,
            activation=nn.ReLU(),
        ):
            """Initialize the neural network.

            Args:
                input_size (int): The size of the input layer.
                hidden_size (int): The size of the hidden layers.
                output_size (int): The size of the output layer.
                num_layers (int): The number of hidden layers.
                dropout (float): The dropout rate.
                activation (nn.Module, optional): The activation function.
                Defaults to nn.ReLU().
            """

            super(NeuralNetModels.Net, self).__init__()
            # make a list of hidden layers
            layers = []
            # input layer
            layers.append(nn.Linear(input_size, hidden_size))
            layers.append(activation)
            layers.append(nn.Dropout(dropout))
            # hidden layers
            for _ in range(1, num_layers):
                layers.append(nn.Linear(hidden_size, hidden_size))
                layers.append(activation)
                layers.append(nn.Dropout(dropout))
            # output layer
            layers.append(nn.Linear(hidden_size, output_size))
            self.net = nn.Sequential(*layers)

        def forward(self, x):
            """Forward pass of the neural network.

            Args:
                x (torch.Tensor): The input data.
            """
            y = self.net(x)
            return y

    def __init__(
        self,
        x_train: np.ndarray,
        y_train: np.ndarray,
        x_test: np.ndarray,
        y_test: np.ndarray,
        save: bool = True,
    ) -> None:
        """Initialize the MaterialsModels class.

        Args:
            x_train (np.ndarray): The training features.
            y_train (np.ndarray): The training labels.
            x_test (np.ndarray): The testing features.
            y_test (np.ndarray): The testing labels.
            save (bool, optional): Whether to save the trained model and
            metrics. Defaults to True.
        """
        super().__init__(x_train, y_train, x_test, y_test, save)

        # convert x and y to torch tensors of type float
        self.x_train = torch.tensor(x_train, dtype=torch.float32)
        self.y_train = torch.tensor(y_train, dtype=torch.float32)
        self.x_test = torch.tensor(x_test, dtype=torch.float32)
        self.y_test = torch.tensor(y_test, dtype=torch.float32)

    def train_models(
        self,
        param_grid: dict = None,
        cv: int = 3,
        seed: int = 42,
        epochs: int = 40,
    ) -> NeuralNetRegressor:
        """Train PyTorch models with cross-validation and grid search.

        Args:
            param_grid (dict): The parameter grid for the grid search.
            Defaults to None.
            cv (int, optional): The number of cross-validation folds.
            Defaults to 3.
            seed (int, optional): The random seed for the train-test split.
            Defaults to 42.
            epochs (int, optional): The number of training epochs.
            Defaults to 40.

        Returns:
            NeuralNetRegressor: The best trained PyTorch model.
        """
        if param_grid is None:
            param_grid = {
                "module__hidden_size": randint(20, 300),
                "module__num_layers": randint(4, 8),
                "module__dropout": uniform(0.0, 0.2),
                "lr": uniform(0.001, 0.01),
            }
        input_size = self.x_train.shape[1]
        output_size = self.y_train.shape[1]

        # set seed for reproducibility
        torch.manual_seed(seed)

        net = NeuralNetRegressor(
            module=self.Net,
            module__input_size=input_size,
            module__hidden_size=15,
            module__output_size=output_size,
            module__num_layers=4,
            module__dropout=0.01,
            max_epochs=epochs,
            lr=0.1,
            optimizer=torch.optim.SGD,
            optimizer__weight_decay=0.0005,
            device="cuda" if torch.cuda.is_available() else "cpu",
            batch_size=256,
        )

        kfold = KFold(n_splits=cv, random_state=seed, shuffle=True)
        random_search = RandomizedSearchCV(
            net,
            param_distributions=param_grid,
            n_iter=20,
            scoring="neg_mean_squared_error",
            n_jobs=3,
            cv=kfold.split(self.x_train, self.y_train),
            verbose=3,
            random_state=seed,
            refit=False,
        )
        random_search.fit(self.x_train, self.y_train)

        # find the best model from the search manually as refit is False
        results = pd.DataFrame(random_search.cv_results_)
        best_idx = results["rank_test_score"].idxmin()
        best_params = results.loc[best_idx, "params"]
        print(f"Best Parameters:\n{best_params}")
        bestnet = NeuralNetRegressor(
            module=self.Net,
            module__input_size=input_size,
            module__hidden_size=best_params["module__hidden_size"],
            module__output_size=output_size,
            module__num_layers=best_params["module__num_layers"],
            module__dropout=best_params["module__dropout"],
            max_epochs=epochs,
            lr=best_params["lr"],
            optimizer=torch.optim.SGD,
            optimizer__weight_decay=0.0005,
            device="cuda" if torch.cuda.is_available() else "cpu",
            batch_size=256,
        )
        bestnet.fit(self.x_train, self.y_train)
        self.model = bestnet

        y_pred = self.model.predict(self.x_test)
        mse = mean_squared_error(self.y_test, y_pred)
        print(f"Mean Squared Error: {mse}")

        if self.save:
            self._dir_out.mkdir(exist_ok=True)
            filename = self._dir_out / "neural_network_model.pkl"
            with open(filename, "wb") as f:
                torch.save(self.model, f)
            pd.DataFrame(best_params, index=[0]).to_csv(
                self._dir_out / "neural_network_params.csv", index=False
            )

        return self.model

    def evaluate_model(self) -> dict:
        """Evaluate the trained model.

        Returns:
            dict: The evaluation metrics.
        """
        super().evaluate_model()

        if self.save:
            self._dir_out.mkdir(exist_ok=True)
            pd.DataFrame(self.metrics, index=[0]).to_csv(
                self._dir_out / "neural_network_metrics.csv", index=False
            )

        return self.metrics
