import sys
import numpy as np
import pandas as pd

import SiGe_model

###****************************** Features and Fingerprints selection *******************
### feature 'a/A' is highly correlating with 'b/A' --> removing 'b/A' (corr.coeff=0.85)
### feature 'alpha' is highly correlating with 'beta' --> removing 'beta' (corr.coeff=0.93)

# possible fingerprints=[
#     'Band_Gap/eV', 
#     'Bulk_modulus',
#     'Thermal_expan_coeff', 
#     'Melting_point',
#     'Thermal_conductivity',
#     'Dielectric_const',
#     'Optical_photon_energy',
#     'Density',
#     'Surface_microhardness',
# ]

features=['a/A', 'c/A', 'alpha', 'gamma', 'volume/A3', 'Si_Num', 'Ge_Num']
not_useful=['formula', 'b/A', 'beta']
###************************** Function calculates number of Si and Ge elements ****
def get_Si_Ge_numbers (df):
    
    formula=df['formula'].str.split('G', n=1, expand=True).rename(columns={0:'Si_Num',1:'Ge_Num'}).fillna(0)
    formula['Si_Num']=formula['Si_Num'].apply(lambda x:0 if x=='' else x)
    formula['Si_Num']=formula['Si_Num'].apply(lambda x:1 if x=='Si' else x)
    formula['Ge_Num']=formula['Ge_Num'].apply(lambda x:1 if x=='e' else x)
    formula['Si_Num']=formula['Si_Num'].apply(lambda x:int(x.split('i')[1]) if type(x)==str else x)
    formula['Ge_Num']=formula['Ge_Num'].apply(lambda x:int(x.split('e')[1]) if type(x)==str else x)
    df['Si_Num']=formula['Si_Num']
    df['Ge_Num']=formula['Ge_Num']
    
    return df
###****************************** Fingerprints Selection for target variables ******************** 
def fprints_selection(train_df):

    not_found=[]
    fprints=[]

    for item in features:
        if item not in list(train_df.columns):
            not_found.append(item)

    if len(not_found)>=1:
        print(f'dataset missing some training features\nPlease add the following for modeling {[i for i in not_found]}, and try again')
        sys.exit()

    for item in list(train_df.columns):
        if (item not in features) & (item not in not_useful):
            fprints.append(item)

    return fprints

###****************************** Make Prediction **********************************    
def estimator(train_df, test_df, fprints):

    for target in fprints:

        y_test_list=[]
        predictions=[]
        R2_list=[]
        model_scores=[]

        trained_model, model_score=SiGe_model.train_model(train_df, target)

        x_test=test_df.loc[:,features]
        pred = trained_model.predict(x_test)
        pred = [round(i,3) for i in pred]
        predictions.append(pred)

        if (target in list(test_df.columns)) & (len(test_df)>2):
            y_test=test_df[target]
            R2=round(trained_model.score(x_test, y_test), 4)
        else:
            R2=model_score
            y_test=[]
        
        y_test_list.append(y_test.to_list())
        R2_list.append(R2)
        model_scores.append(model_score)

    return y_test_list, predictions, R2_list, model_scores

###****************************** Tune the model **********************************   
def tune(train_df, test_df, fprints):
    
    for target in fprints:
        y_test_list=[]
        predictions=[]
        R2_list=[]
        model_scores=[]
        
        tuned_model, tuned_model_score=SiGe_model.tune_model(train_df, test_df, target)

        x_test=test_df.loc[:,features]
        pred = tuned_model.predict(x_test)
        pred = [round(i,3) for i in pred]
        predictions.append(pred)

        if (target in list(test_df.columns)) & (len(test_df)>2):
            y_test=test_df[target]
            R2=round(tuned_model.score(x_test, y_test), 4)
        else:
            R2=tuned_model_score
            y_test=[]

        y_test_list.append(y_test.to_list())
        R2_list.append(R2)
        model_scores.append(tuned_model_score)

    return y_test_list, predictions, R2_list, model_scores


