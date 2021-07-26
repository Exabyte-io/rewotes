###****************************** Model_v2 ***************************************
### Simple version for now
###****************************** Libraries Imports ******************************
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.svm import SVR 
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from scipy.stats import randint

###****************************** Features and Finger Prints******************************
features=['a/A', 'b/A', 'c/A', 'alpha', 'beta', 'gamma', 'volume/A3', 'Si_Num', 'Ge_Num']
finger_prints=['Band_Gap/eV']

###****************************** Model Train (Pipleline)**********************************
def train_model(df):
#   print('data size: ', df.shape)

    X=df.loc[:,features]
    y=df[finger_prints[0]]

    X_train, X_test, y_train, y_test=train_test_split(X, y, test_size=0.1, random_state=775)

    pipe=Pipeline(steps=[('mmscaler', MinMaxScaler()), ('svr', SVR(C=2.2, epsilon=0.18))])
    pipe.fit(X_train, y_train)  #tuning is planned with Grid/Random Search later on
    trained_model = pipe
    model_score=round(pipe.score(X_test, y_test),4)

    return trained_model, model_score

###****************************** Model Tune (Randomized Search CV)**********************************
def tune_model(df, test_df):

    X=df.loc[:,features]
    y=df[finger_prints[0]]

    if finger_prints[0] in list(test_df.columns):
        y_test=test_df[finger_prints[0]]
        y_train=y
        X_test=test_df.loc[:,features]
        X_train=X
    else:
        X_train, X_test, y_train, y_test=train_test_split(X, y, test_size=0.2, random_state=573)

    C_list=[float(x) for x in np.linspace(3.0, 1.5, num=24)]
    epsilon_list=[float(x) for x in np.linspace(0.15, 0.35, num=24)]

    params={
        'svr__kernel':['rbf', 'linear', 'poly', 'sigmoid'],
        'svr__gamma':['scale', 'auto'],
        'svr__degree':randint(2,8),
        'svr__C':C_list,
        'svr__epsilon':epsilon_list
    }
    pipe=Pipeline(steps=[('mmscaler', MinMaxScaler()), ('svr', SVR())])

    samples=150
    randomCV=RandomizedSearchCV(pipe, param_distributions=params, n_iter=samples)
    randomCV.fit(X_train, y_train)  
    tuned_model = randomCV
    tuned_model_score=round(randomCV.score(X_test, y_test),4)

    return tuned_model, tuned_model_score
    


