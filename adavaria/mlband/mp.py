import pandas as pd
import warnings
from pymatgen.io.cif import CifWriter


def get_mp_api_connector():
    from mp_api.client import MPRester
    api_key = 'pJmFwIk6qWCnMwEuGvzwLBT1qfdPpElH'
    return MPRester(api_key)

mpr = None
try:    
    mpr = get_mp_api_connector()
except:
    warnings.warn('Could not connect to MP API.')
    pass

def mp_summary_to_df(data, drop_fields_not_requested=True):
    df = pd.DataFrame([doc.dict() for doc in data])
    if drop_fields_not_requested:
        df.drop(columns=['fields_not_requested'], inplace=True)
    return df