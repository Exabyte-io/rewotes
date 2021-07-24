import os
import numpy as np
import pandas as pd

from matminer.utils.io import load_dataframe_from_json

from sklearn.model_selection import train_test_split

from megnet.data.crystal import CrystalGraph
from megnet.data.graph import GaussianDistance
from megnet.models import MEGNetModel

from utils import clean_dataframe

DATA_SEED = 42
R_CUT = 4
R_BASIS = 10
R_SMEAR = 0.5
PATH = os.path.dirname(os.path.abspath(__file__))

data = "non-metals"
# data = "combined"

df = load_dataframe_from_json(
    # PATH + "/data/mp-query-exabyte_24_07_2021.json.gz"
    PATH + "/data/mp-query-exabyte-test_23_07_2021.json.gz"
)

print("Clean the Data Set")
df = clean_dataframe(df, R_CUT)


print("Perform Stratified Train-Test Split")
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


converter = CrystalGraph(bond_converter=GaussianDistance(np.linspace(0, R_CUT, R_BASIS), R_SMEAR))

model = MEGNetModel(R_BASIS,
    graph_converter=converter
)

model.train(df_train["final_structure"], df_train["band_gap"], epochs=100)

model.save_model(PATH+f"/models/structure-megnet-{data}.hdf5")
