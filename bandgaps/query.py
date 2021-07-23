import os
import time
import pandas as pd

from pymatgen.core import SETTINGS
from pymatgen.ext.matproj import MPRester
from matminer.utils.io import store_dataframe_as_json

PATH = os.path.dirname(os.path.abspath(__file__))
# Materials project provides a python API via pymatgen
# API key obtained from https://materialsproject.org/open
# place the key in ~/.pmgrc.yaml config file
MP_API_KEY = SETTINGS.get("PMG_MAPI_KEY")

m = MPRester(MP_API_KEY)

# NOTE magmoms cannot be queried in this way
target_properties = [
    "band_gap",
    "e_above_hull",
    "final_structure",
]

# Only take data where the bandstructure has been calculated directly
# this avoids issues due to band_gaps calculated with loose k-grids
# target_criteria = {"has": "bandstructure"}  # 75067 in total
target_criteria = {
    "e_above_hull": {"$lt": 0.01},
    "has": "bandstructure",
    "nelements": 2,
}  # 5909 use for development

data = m.query(criteria=target_criteria, properties=target_properties)

df = pd.DataFrame(data)

# Save dataframe storing MSONable objects as json i.e. pymatgen Structures
date = time.strftime("%d_%m_%Y")
store_dataframe_as_json(
    df, PATH + f"/data/mp-query-exabyte-test_{date}.json", compression="gz"
)
