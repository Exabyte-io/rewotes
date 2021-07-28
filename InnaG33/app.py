###*********************** Performs estimations of SiGe properties based on the lattice parameters *********
'''
> Minimum required features are SiGe formula (or Si, Ge alone), lattice size, angles, and volume
> User's train and test dataset are used from 'data' subfolder
> Automated selection of minimum required features from user's train / test dataset
> Automated alocation for rest of found features as target varibales (SiGe fingerprints)
> SVR algorithm deals with small datasets (suggested data size up to 10,000 rows)
> Option for model tune --> per user's choice and data quality
> Test dataset (with min one row) can only include minimum required features --> 
model will predict target variable per what is found in the train dataset
and model score will be based on a split from the train dataset.
'''

import sys
import pandas as pd

import SiGe_estimator

###****************************** Data read ***************************************
try:
    train_df=pd.read_csv('./data/SiGe_lattice_params_train.csv', low_memory=False, encoding = 'utf8') 
    test_df=pd.read_csv('./data/SiGe_lattice_params_test.csv',low_memory=False, encoding = 'utf8')   #separate test data file (avoid data leakage on model train)
except: 
    FileNotFoundError
    print('please add your dataset in the "data" folder and try again')
    sys.exit()

###****************************** Data prep for modeling **************************
train_df=SiGe_estimator.get_Si_Ge_numbers(train_df)
train_df.drop(['formula'], axis=1, inplace=True)

test_df=SiGe_estimator.get_Si_Ge_numbers(test_df)
test_df.drop(['formula'], axis=1, inplace=True)

###****************************** Fingerprints Selection for target variables *************
fprints=SiGe_estimator.fprints_selection(train_df)
print()
print('Found following fprints for modeling: ', fprints)
print()
###****************************** Train the model and make prediction *************
y_test_list, predictions, R2_list, model_scores = SiGe_estimator.estimator(train_df, test_df, fprints)

for i in range(len(fprints)):
    print(
        f'Model score is: {model_scores[i]}\nActual {fprints[i]}: {y_test_list[i]}\nPredicted {fprints[i]}: {predictions[i]}\nTest data score {R2_list[i]}'
        )
    print()
###****************************** Model tune - optional per user's data & choice *************

user_choice=str(input('Would you like try tuning the model (may take few minutes): y/n? '))
print()

if user_choice.lower()=='y':
    y_test_list, predictions, R2_list, model_scores = SiGe_estimator.tune(train_df, test_df, fprints)
    for i in range(len(fprints)):
        print(
            f'Model score is: {model_scores[i]}\nActual {fprints[i]}: {y_test_list[i]}\nPredicted {fprints[i]}: {predictions[i]}\nTest data score {R2_list[i]}'
            )
        print()
else:
    print ('run completed')

print()