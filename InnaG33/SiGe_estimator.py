import numpy as np
import pandas as pd
import random

import SiGe_model

###****************************** Features and Fingerprints selection *******************
### feature 'a/A' is highly correlating with 'b/A' --> removing 'a/A' (corr.coeff=0.85)
### feature 'alpha' is highly correlating with 'beta' --> removing 'alpha' (corr.coeff=0.93)

features=['b/A', 'c/A', 'beta', 'gamma', 'volume/A3', 'Si_Num', 'Ge_Num']
fingerprints=[
    'Band_Gap/eV', 
    'Bulk_modulus',
    'Thermal_expan_coeff', 
    'Melting_point',
    'Thermal_conductivity',
    'Dielectric_const',
    'Optical_photon_energy',
    'Density',
    'Surface_microhardness',
    'Number_atoms_1cm3'
]

###****************************** Make Prediction **********************************    
def estimator(train_df, test_df):

    for i in range(len(fingerprints)):
        if fingerprints[i] in list(train_df.columns):
            target=fingerprints[i]
#            updated_fp=fingerprints[i+1:]
            break

    trained_model, model_score=SiGe_model.train_model(train_df, target)

    x_test=test_df.loc[:,features]
    pred = trained_model.predict(x_test)
    pred = [round(i,3) for i in pred]

    if target in list(test_df.columns):
        y_test=test_df[target]
        R2=round(trained_model.score(x_test, y_test), 4)
    else:
        R2=model_score
        y_test=[]

    return y_test, pred, R2, model_score

def tune(train_df, test_df):
    
    for i in range(len(fingerprints)):
        if fingerprints[i] in list(train_df.columns):
            target=fingerprints[i]
#            updated_fp=fingerprints[i+1:]
            break
    
    tuned_model, tuned_model_score=SiGe_model.tune_model(train_df, test_df, target)

    x_test=test_df.loc[:,features]
    pred = tuned_model.predict(x_test)
    pred = [round(i,3) for i in pred]

    if target in list(test_df.columns):
        y_test=test_df[target]
        R2=round(tuned_model.score(x_test, y_test), 4)
    else:
        R2=tuned_model_score
        y_test=[]

    return y_test, pred, R2, tuned_model_score


