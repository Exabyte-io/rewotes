# Train a composition-based model to predict bandgaps
import os
import pandas as pd
from joblib import dump

from matminer.featurizers import composition as cf
from matminer.featurizers.base import MultipleFeaturizer
from matminer.utils.io import load_dataframe_from_json

from sklearn.model_selection import train_test_split

from utils import RandomForestRegressor


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
dump(model, PATH+f"/models/composition-rfr-{data}.pkl.gz", )
