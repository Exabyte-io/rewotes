import os
import numpy as np
import pandas as pd
from joblib import load

import matplotlib.pyplot as plt

from pymatgen.core.composition import Composition

from matminer.featurizers import composition as cf
from matminer.featurizers.base import MultipleFeaturizer

from utils import RandomForestRegressor  # needed to load model


PATH = os.path.dirname(os.path.abspath(__file__))


# data = "non-metals"
data = "combined"

model = load(PATH+f"/models/composition-rfr-{data}.pkl.gz")

comp_features = MultipleFeaturizer([
    cf.Stoichiometry(),
    cf.ElementProperty.from_preset("magpie"),
    cf.ValenceOrbital(props=["avg"]),
    cf.IonProperty(fast=True)
])

feature_labels = comp_features.feature_labels()

print("Evaluate results on Si-Ge chemical System")
df_test_comps = pd.DataFrame()
df_test_comps["composition"] = [Composition({"Si": 1-x, "Ge": x}) for x in np.linspace(0, 1)]
df_test_comps = comp_features.featurize_dataframe(df_test_comps, col_id="composition")

# y_pred = model.predict(df_test_comps[feature_labels])
y_pred, y_unc = model.predict(df_test_comps[feature_labels], return_std=True)

plt.figure(figsize=(8, 8))
plt.xlim
plt.plot(np.linspace(0, 1), y_pred, label="RFR Prediction")
plt.fill_between(np.linspace(0, 1), y_pred+y_unc, y_pred-y_unc, alpha=0.4)

# https://journals.aps.org/prb/abstract/10.1103/PhysRevB.40.5683
fitted = [1.155 - 0.43*x + 0.206*x**2 if x < 0.85 else 2.010 - 1.270*x for x in np.linspace(0, 1)]

plt.plot(np.linspace(0, 1), fitted, label="Experimental Fit")
ax = plt.gca()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.5)

# plt.scatter((0, 0.125, 0.25, 0.5, 0.75, 1), (0.8527, 0.8585, 0.0000, 0.4267, 0.0000, 0.0000), color="red", marker="x", clip_on=False, zorder=3, label="Training Data")
# plt.scatter((0, 0.125, 0.5), (0.8527, 0.8585, 0.4267), color="red", marker="x", clip_on=False, zorder=3, label="Training Data")

plt.legend(frameon=False)

xlabel = r"$x$" + " in " + r"$Si_{1-x}Ge_{x}$"
ax.set_xlabel(xlabel)

ylabel = r"$E_g$" + " / eV"
ax.set_ylabel(ylabel)
plt.savefig(PATH+f"/SiGe-composition-{data}.pdf")
plt.show()
