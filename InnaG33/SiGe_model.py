###****************************** Model_v0 ***************************************
### Simplest version, planned to separate into "model" and "estimator" run by one "app.py"

###****************************** Libraries Imports ******************************
import numpy as np
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.svm import SVR 

###****************************** Data read ***************************************
df=pd.read_csv('SiGe_lattice_params_BG.csv')

X=df.iloc[:,1:10]   #features selection is planned
y=df['Band_Gap/eV'] #finger prints selection is planned

###****************************** Model Pipeline --> Simplest model for v0 ********
### No train-test-valid split due very small data size

pipe=Pipeline(steps=[('mmscaler', MinMaxScaler()), ('svr', SVR(C=2.0, epsilon=0.02))])
pipe.fit(X, y)  #tuning is planned with GridSearch later on
print('model score: ', round(pipe.score(X, y),3))

user_test = int(input('eneter row number to test, from 0 to 5: ')) #preliminary, test can be loaded from separate csv

x_test=df.iloc[user_test,1:10]
y_test=df.iloc[user_test,10]

pred_1 = pipe.predict(x_test.values.reshape(1,-1))
print(f'Actual Bang Gap value: {y_test}, Predicted: {round(pred_1[0],3)}\nResidual is {round(np.sqrt((y_test-pred_1[0])**2), 4)}')


