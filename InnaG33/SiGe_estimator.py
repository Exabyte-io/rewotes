import numpy as np
import pandas as pd

import SiGe_model

###****************************** Make Prediction **********************************    
def estimator(train_df, test_df):

    features=['a/A', 'b/A', 'c/A', 'alpha', 'beta', 'gamma', 'volume/A3', 'Si_Num', 'Ge_Num']
    finger_prints=['Band_Gap/eV']

    x_test=test_df.loc[:,features]
    y_test=test_df[finger_prints[0]]  #if y_test provided, else: no test score would be avaialble (remove R2 from model output) 

    trained_model, model_score=SiGe_model.train_model(train_df)

    pred = trained_model.predict(x_test)
    pred = [round(i,3) for i in pred]

    R2=round(trained_model.score(x_test, y_test), 4)

    return y_test, pred, R2, model_score
