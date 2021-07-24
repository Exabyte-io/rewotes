# Train a composition-based model to predict bandgaps
import os
import numpy as np
import pandas as pd
from joblib import dump

import matplotlib.pyplot as plt

from pymatgen.core.composition import Composition

from matminer.featurizers import composition as cf
from matminer.featurizers.base import MultipleFeaturizer
from matminer.utils.io import load_dataframe_from_json

from sklearn.model_selection import train_test_split
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


DATA_SEED = 42
R_CUT = 4
R_BASIS = 10
R_SMEAR = 0.5
PATH = os.path.dirname(os.path.abspath(__file__))

data = "non-metals"
# data = "combined"

df = load_dataframe_from_json(
    PATH + "/data/mp-query-exabyte_24_07_2021.json.gz"
    # PATH + "/data/mp-query-exabyte-test_23_07_2021.json.gz"
)

print("Clean the Data Set")
# As we're using a composition based model here we cannot distinguish polymorphs
# to ensure a clear mapping between composition and the target
df["composition"] = df["final_structure"].apply(lambda x: x.composition)
df["formula"] = df["composition"].apply(lambda x: x.reduced_formula)
df = df.sort_values(by=["formula", "e_above_hull"], ascending=True, kind="mergesort")
df = df.drop_duplicates(["formula"], keep="first")

# discard structures with isolated atoms (no neighbours within 4\AA)
all_iso = []
some_iso = []
for idx, crystal in zip(df.index, df["final_structure"]):
    self_idx, nbr_idx, *_ = crystal.get_neighbor_list(
        R_CUT,
        numerical_tol=1e-8,
    )

    if len(self_idx) == 0:
        all_iso.append(idx)
    elif len(nbr_idx) == 0:
        all_iso.append(idx)
    elif set(self_idx) != set(range(crystal.num_sites)):
        some_iso.append(idx)

if (len(all_iso) > 0) or (len(some_iso) > 0):
    # drop the structures with isolated atoms
    df = df.drop(df[df.index.isin(all_iso)].index)
    df = df.drop(df[df.index.isin(some_iso)].index)

# Use the features from MAGPIE
print("Construct Composition-Based Descriptors")
comp_features = MultipleFeaturizer([
    cf.Stoichiometry(),
    cf.ElementProperty.from_preset("magpie"),
    cf.ValenceOrbital(props=["avg"]),
    cf.IonProperty(fast=True)
])

feature_labels = comp_features.feature_labels()
df = comp_features.featurize_dataframe(df, col_id="composition")

df = df.dropna()

print("Perform Stratefied Train-Test Split")
df_metals = df[df["band_gap"] == 0]
df_non_metals = df[df["band_gap"] > 0]

df_metal_train, df_metal_test = train_test_split(df_metals, test_size=0.2, random_state=DATA_SEED)
df_non_metal_train, df_non_metal_test = train_test_split(df_non_metals, test_size=0.2, random_state=DATA_SEED)

if data == "combined":
    df_train = pd.concat((df_metal_train, df_non_metal_train))
    df_test = pd.concat((df_metal_test, df_non_metal_test))
elif data == "non-metals":
    df_train = df_non_metal_train
    df_test = df_non_metal_test

model = RandomForestRegressor()

print("Train Random Forest Regressor")
model.fit(df_train[feature_labels].values, df_train["band_gap"].values)

# Save the model
dump(model, PATH+f"/models/composition-rfr-{data}.pkl")

print("Evaluate results on Si-Ge chemical System")
df_test_comps = pd.DataFrame()
df_test_comps["composition"] = [Composition({"Si": 1-x, "Ge": x}) for x in np.linspace(0, 1)]
df_test_comps = comp_features.featurize_dataframe(df_test_comps, col_id="composition")

# y_pred = model.predict(df_test_comps[feature_labels])
y_pred, y_unc = model.predict(df_test_comps[feature_labels], return_std=True)

plt.figure(figsize=(8, 8))
plt.xlim
plt.plot(np.linspace(0, 1), y_pred, label="RFR Prediction")
plt.fill_between(np.linspace(0, 1), y_pred+y_unc, y_pred-y_unc, label="RFR Prediction", alpha=0.4)

# https://journals.aps.org/prb/abstract/10.1103/PhysRevB.40.5683
fitted = [1.155 - 0.43*x + 0.206*x**2 if x < 0.85 else 2.010 - 1.270*x for x in np.linspace(0, 1)]

plt.plot(np.linspace(0, 1), fitted, label="Experimental Fit")
ax = plt.gca()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.5)

# plt.scatter((0, 0.125, 0.25, 0.5, 0.75, 1), (0.8527, 0.8585, 0.0000, 0.4267, 0.0000, 0.0000), color="red", marker="x", clip_on=False, zorder=3, label="Training Data")
# plt.scatter((0, 0.125, 0.5), (0.8527, 0.8585, 0.4267), color="red", marker="x", clip_on=False, zorder=3, label="Training Data")

plt.legend(frameon=False)

xlabel = r"$x$" + " in " + r"$Si_{x}Ge_{1-x}$"
ax.set_xlabel(xlabel)

ylabel = r"$E_g$" + " / eV"
ax.set_ylabel(ylabel)
plt.savefig(PATH+f"/SiGe-composition-{data}.pdf")
plt.show()
