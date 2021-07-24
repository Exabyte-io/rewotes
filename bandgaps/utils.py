import numpy as np

from sklearn.ensemble import RandomForestRegressor as skrfr

class RandomForestRegressor(skrfr):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def predict(self, X_test, return_std=False):
        """Predict continuous labels y_pred and uncertainty y_std
        (unless no_std=True) for X_test.

        Args:
            X_test (array-like, shape=(n_samples, n_features)): Input data.
            no_std=False (bool, optional): Don't return y_std if set to true.

        Returns:
            2- or 1-tuple: y_pred (and y_std)
        """
        if self.criterion != "mse":
            err = f"Expected impurity to be 'mse', instead got {self.criterion}."
            raise ValueError(err)

        y_pred = super().predict(X_test)
        if return_std:
            y_std = self._return_std(X_test, self.estimators_, y_pred)
            return y_pred, y_std
        else:
            return y_pred


    @staticmethod
    def _return_std(X_test, trees, predictions):
        """Returns `std(Y | X_test)` calculated via E[Var(Y | Tree)] + Var(E[Y | Tree])
        where P(Tree) is `1 / len(trees)`.

        Note: Another option for estimating confidence intervals would be prediction
        variability, i.e. how influential training set is for producing observed random
        forest predictions. E.g. https://github.com/scikit-learn-contrib/forest-confidence-interval.
        Empirically the below method to obtain y_std was found to be noticeably more accurate.

        Args:
            X_test (array-like, shape=(n_samples, n_features)): Input data.
            trees (list, shape=(n_estimators,)): List of fitted sklearn trees as obtained
                from the estimators_ attribute of a fitted RandomForestRegressor.
            predictions (array-like, shape=(n_samples,)): Prediction for each sample
                as returned by RandomForestRegressor.

        Returns:
            array-like, shape=(n_samples,): Standard deviation of y_pred at X_test. If criterion
            is set to "mse", then std[i] ~= std(y | X_test[i]).
        """
        y_var = np.zeros(len(X_test))

        # This derives std(y | X_test) as described in sec. 4.3.2 of
        # http://arxiv.org/pdf/1211.0906v2.pdf
        for tree in trees:
            var_tree = tree.tree_.impurity[tree.apply(X_test)]
            mean_tree = tree.predict(X_test)
            y_var += var_tree + mean_tree ** 2

        y_var /= len(trees)
        y_var -= predictions ** 2
        return y_var ** 0.5