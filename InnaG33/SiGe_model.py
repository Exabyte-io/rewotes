###****************************** Model_v1 ***************************************
### Simplest version, planned to separate into "model" and "estimator" run by one "app.py"

###****************************** Libraries Imports ******************************
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.svm import SVR 

#import SiGe_model
#import SiGe_estimator

###****************************** Model Train (Pipleline)**********************************
def train_model(df):

    X=df.iloc[:,1:10]   #features selection is planned
    y=df['Band_Gap/eV'] #finger prints selection is planned

### No train-test-valid split due very small data size, so that score based on train data
    pipe=Pipeline(steps=[('mmscaler', MinMaxScaler()), ('svr', SVR(C=2.2, epsilon=0.06))])
    pipe.fit(X, y)  #tuning is planned with GridSearch later on
    trained_model = pipe
    model_score=round(pipe.score(X, y),4)

    return trained_model, model_score
    
###****************************** Make Prediction **********************************    
def estimator(train_df, test_df):

    x_test=test_df.iloc[:,1:10]
    y_test=test_df['Band_Gap/eV']

    trained_model, model_score=train_model(train_df)

    pred = trained_model.predict(x_test)
    pred = [round(i,3) for i in pred]

    residuals=sum([(y_test[i]-pred[i])**2 for i in range(len(y_test))])
    R2=round(np.sqrt(residuals),4)

    return y_test, pred, R2, model_score

###****************************** Data read ***************************************
train_df=pd.read_csv('SiGe_lattice_params_train.csv') 
test_df=pd.read_csv('SiGe_lattice_params_test.csv')   #separate test data file (avoid data leakage on model train)

###****************************** Train the model and make prediction *************
y_test, pred, R2, model_score = estimator(train_df, test_df)

print(f'Model score is: {model_score}\nActual Bang Gap: {y_test.to_list()}\nPredicted Band Gap: {pred}\nResiduals {R2}')

