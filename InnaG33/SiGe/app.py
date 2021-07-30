###*********************** Performs estimations of SiGe properties based on the lattice parameters *********
'''
> Minimum required features are SiGe formula (or Si, Ge alone), lattice size, angles, and volume
> User's train and test dataset are used from 'data' subfolder
> Automated selection of minimum required features from user's train / test dataset
> Automated allocation for rest of found features as target varibales (SiGe fingerprints)
> SVR algorithm deals with small datasets (suggested data size up to 10,000 rows)
> Option for model tune 
> Test dataset (with min one row) can only include minimum required features --> 
model will predict target variable per what is found in the train dataset
and model score will be based on a split from the train dataset.
'''

from flask import Flask, jsonify
import sys
import pandas as pd
import SiGe_estimator

#################################################
# Flask Setup
#################################################
app = Flask(__name__)

###****************************** Data read ***************************************
try:
    path='/Users/grinn/Desktop/rewotes/InnaG33/SiGe'
    train_df=pd.read_csv(f'{path}/data/train_df.csv', low_memory=False, encoding = 'utf8') 
    test_df=pd.read_csv(f'{path}/data/test_df.csv',low_memory=False, encoding = 'utf8')   #separate test data file (avoid data leakage on model train)
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

#################################################
# Flask Routes
#################################################

@app.route("/")
def welcome():
    """List all available api routes."""
    return (
        f"Available Routes:<br/>"
        f"/model_train<br/>"
        f"/model_tune<br/>"
    )

###****************************** Train the model and make prediction ********************
@app.route("/model_train")
def model_train():

    results={}

    y_test_list, predictions, R2_list, model_scores = SiGe_estimator.estimator(train_df, test_df, fprints)

    for i in range(len(fprints)):

        results[f"Model score for {fprints[i]}"] = model_scores[i]
        results[f"Actual {fprints[i]}"] = y_test_list[i]
        results[f"Predicted {fprints[i]}"] = predictions[i]
        results[f"Test data score {fprints[i]}"] = R2_list[i]

    return jsonify(results)

###****************************** Model tune - optional per user's data & choice *************

@app.route("/model_tune")
def model_tune():
    results={}

    y_test_list, predictions, R2_list, model_scores = SiGe_estimator.tune(train_df, test_df, fprints)

    for i in range(len(fprints)):
        results[f"Model score for {fprints[i]}"] = model_scores[i]
        results[f"Actual {fprints[i]}"] = y_test_list[i]
        results[f"Predicted {fprints[i]}"] = predictions[i]
        results[f"Test data score {fprints[i]}"] = R2_list[i]

    return jsonify(results)

if __name__ == '__main__':
    app.run(debug=True)