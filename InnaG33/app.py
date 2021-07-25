import pandas as pd

import SiGe_estimator

###************************** Function calculates number of Si and Ge elements ****
### can be transformed to Preprocess class****
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

###****************************** Data read ***************************************
train_df=pd.read_csv('SiGe_lattice_params_train.csv') 
test_df=pd.read_csv('SiGe_lattice_params_test.csv')   #separate test data file (avoid data leakage on model train)

###****************************** Data prep for modeling **************************
train_df=get_Si_Ge_numbers(train_df)
train_df.drop(['formula'], axis=1, inplace=True)

test_df=get_Si_Ge_numbers(test_df)
test_df.drop(['formula'], axis=1, inplace=True)

###****************************** Train the model and make prediction *************

y_test, pred, R2, model_score = SiGe_estimator.estimator(train_df, test_df)

print(f'Model score is: {model_score}\nActual Bang Gap: {y_test.to_list()}\nPredicted Band Gap: {pred}\nTest data score {R2}')