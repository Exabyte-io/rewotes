###****************************** Model_v2 ***************************************
### Simple version for now
###****************************** Libraries Imports ******************************
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.svm import SVR 
from sklearn.model_selection import train_test_split

###****************************** Model Train (Pipleline)**********************************
def train_model(df):

    features=['a/A', 'b/A', 'c/A', 'alpha', 'beta', 'gamma', 'volume/A3', 'Si_Num', 'Ge_Num']
    finger_prints=['Band_Gap/eV']
#    print('data size: ', df.shape)

    X=df.loc[:,features]
    y=df[finger_prints[0]]

    X_train, X_test, y_train, y_test=train_test_split(X, y, test_size=0.2, random_state=773)

    pipe=Pipeline(steps=[('mmscaler', MinMaxScaler()), ('svr', SVR(C=2.2, epsilon=0.18))])
    pipe.fit(X_train, y_train)  #tuning is planned with Grid/Random Search later on
    trained_model = pipe
    model_score=round(pipe.score(X_test, y_test),4)

    return trained_model, model_score
    


