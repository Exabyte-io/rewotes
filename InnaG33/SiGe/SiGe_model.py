###****************************** Model_v4 with Features Fingerprints automated selection **************
### SVR for smaller datasets
###****************************** Libraries Imports ******************************
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.svm import SVR 
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from scipy.stats import randint

###****************************** Features and Fingerprints selection *******************
### feature 'a/A' is highly correlating with 'b/A' --> removing 'a/A' (corr.coeff=0.85)
### feature 'alpha' is highly correlating with 'beta' --> removing 'alpha' (corr.coeff=0.93)

features=['b/A', 'c/A', 'beta', 'gamma', 'volume/A3', 'Si_Num', 'Ge_Num']

###****************************** Model Train (Pipleline)**********************************
def train_model(df, target):

    X=df.loc[:,features]
    y=df[target]

    X_train, X_test, y_train, y_test=train_test_split(X, y, test_size=0.1, random_state=775)

    pipe=Pipeline(steps=[('mmscaler', MinMaxScaler()), ('svr', SVR(C=2.2, epsilon=0.18))])
    pipe.fit(X_train, y_train)  
    trained_model = pipe
    model_score=round(pipe.score(X_test, y_test),4)

    return trained_model, model_score

###****************************** Model Tune (Randomized Search CV)**********************************
def tune_model(df, test_df, target):

    X=df.loc[:,features]
    y=df[target]

    if target in list(test_df.columns):
        y_test=test_df[target]
        y_train=y
        X_test=test_df.loc[:,features]
        X_train=X
    else:
        X_train, X_test, y_train, y_test=train_test_split(X, y, test_size=0.1, random_state=573)

    C_list=[float(x) for x in np.linspace(1.5, 4.0, num=12)]
    epsilon_list=[float(x) for x in np.linspace(0.15, 0.4, num=12)]

    params={
        'svr__kernel':['rbf', 'linear', 'poly', 'sigmoid'],
        'svr__gamma':['scale', 'auto'],
        'svr__degree':randint(2,6),
        'svr__C':C_list,
        'svr__epsilon':epsilon_list
    }
    pipe=Pipeline(steps=[('mmscaler', MinMaxScaler()), ('svr', SVR())])

    samples=200
    randomCV=RandomizedSearchCV(pipe, param_distributions=params, n_iter=samples, cv=3)
    randomCV.fit(X_train, y_train)  
    tuned_model = randomCV
    tuned_model_score=round(randomCV.score(X_test, y_test),4)

    return tuned_model, tuned_model_score
    


