import numpy as np
import pandas as pd
import random

import SiGe_model

features=['a/A', 'b/A', 'c/A', 'alpha', 'beta', 'gamma', 'volume/A3', 'Si_Num', 'Ge_Num']
finger_prints=['Band_Gap/eV']

###****************************** Make Prediction **********************************    
def estimator(train_df, test_df):

    trained_model, model_score=SiGe_model.train_model(train_df)

    x_test=test_df.loc[:,features]
    pred = trained_model.predict(x_test)
    pred = [round(i,3) for i in pred]

    if finger_prints[0] in list(test_df.columns):
        y_test=test_df[finger_prints[0]]
        R2=round(trained_model.score(x_test, y_test), 4)
    else:
        R2=model_score

    return y_test, pred, R2, model_score

def tune(train_df, test_df):

    tuned_model, tuned_model_score=SiGe_model.tune_model(train_df, test_df)

    x_test=test_df.loc[:,features]
    pred = tuned_model.predict(x_test)
    pred = [round(i,3) for i in pred]

    if finger_prints[0] in list(test_df.columns):
        y_test=test_df[finger_prints[0]]
        R2=round(tuned_model.score(x_test, y_test), 4)
    else:
        R2=tuned_model_score

    return y_test, pred, R2, tuned_model_score


