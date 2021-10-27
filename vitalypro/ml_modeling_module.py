# -*- coding: utf-8 -*-
"""
@author: Vitaly Proshchenko
"""

from datetime import datetime        
from hyperopt import hp, fmin, tpe, STATUS_OK, Trials, space_eval
from hyperopt.fmin import generate_trials_to_calculate
from sklearn.model_selection import cross_val_score                
from functools import partial
import numpy as np
        
   
    
"""
1. Add class for data (energy bands, cell geometry etc) extraction from QE output files
2. Add class for data preparation for ML modeling
3. Add class for ML modeling with RF, XGB (more ML algorithms -?) 
4. If time permitted, add class to solve reverse problems (e.g. for known bandstructure, find SL structure)
""" 

class modeling:
    def __init__(self, X_train, Y_train, X_test, Y_test, scaler= 'None'):
        """
        since we use tree based ML algorithms,
        there is no need to scale data (scaler='None').
        In the future work scaler='minmax', 'standard' etc can be added for other ML algorithms
          
        """
        self.scaler = scaler
        if scaler== 'None':
            self.X_train = X_train
            self.Y_train = Y_train
            self.X_test = X_test
            self.Y_test = Y_test
    
        """ if scaler != 'None' -- perform data scaling"""  
    
     
    def verbose_output(self,  model,  y_fit, y_predict):
        """ print train/test model scores if verbose output is chosen"""
        
        from sklearn.metrics import mean_squared_error
        from sklearn.metrics import r2_score

        text=("\nTrain:\nR^2:"+str(round(model.score(self.X_train, self.Y_train),2))+
            ' RMSE:'+str(round(np.sqrt(mean_squared_error(self.Y_train, y_fit, multioutput='uniform_average')),2))+
            '\nTest:\nR^2:'+str(round(r2_score(self.Y_test, y_predict, multioutput='uniform_average'),2)) + 
            ' RMSE: '+str( round(np.sqrt(mean_squared_error(self.Y_test, y_predict)),2))  
            )
        return text
    
    def RandomForest(self,  max_evals=40, timeout=600,  loss_threshold=0.5, init_vals='None' , n_jobs=10, cv=5, verbose=1):
        """ 
        Random Forest Regressor with Sequential Model-based Optimization with Tree-Parzen estimators (Bayesian optimization) 
        The idea is to build a probability model of the objective function and use it to select the most promising hyperparameters to evaluate in the true objective function
        
        :max_evals - number of max evaluation in hyperparameter search algorithm. the search will stop when number of max itteration  is reached. Can be any integer number. Default value is 40, 
        :timeout - search time. the search will stop when the time is out. Can be any integer number. Default value is 600 sec  
        :loss_threshold - the search will stop when the loss threshold is achieved. (in our case loss function is MSE). Can be any float number.
        :init_vals - initial guess for the hyperparameter optimization search. Default - 'None'. Can be list of parameters such as n_estimators, max_depth etc.  
        :n_jobs - number of CPUs for parallelization = 10, 
        :cv - cross-validation. default = 5 
        :verbose - verbose output. Can be either 0 or 1.
        
        """
        from sklearn.ensemble import RandomForestRegressor        
        startTime = datetime.now()
    
        def objective(space):
            """Define our function"""
            clf = RandomForestRegressor(n_estimators =int(space['n_estimators']),                           
                                        random_state = int(space['random_state']),
                                        max_depth = int(space['max_depth']),
                                        min_samples_split = int(space['min_samples_split'])
                                        )
    
            """To reduce overfitting, we perform n-fold cross validation.
                The goal is to minimize our loss function (MSE)"""
            mse_scr = cross_val_score(clf, self.X_train, self.Y_train, scoring='neg_root_mean_squared_error', cv=cv,verbose=0,n_jobs=n_jobs).mean()
    
            return {'loss': -mse_scr, 'status': STATUS_OK }
        
        """Here we define the space for Random Forest model
           More hyperparameters to build the sapce can also be used """    
        
        space ={'max_depth':  hp.quniform('max_depth', 1, 18, 1),          
                'random_state' : hp.quniform('random_state', 1, 18, 1),
                'n_estimators' : hp.quniform('n_estimators', 10, 800, 1),
                'min_samples_split' : hp.quniform('min_samples_split', 1, 8, 1)
                }
        
        """ in case if we do not provide the initial guess for Random forest model, 
            we will first randomly pick 8 sets of hyperparameters (HPs) and build a probabalistic space. 
            The next HPs will be chosen based on probabalistic model
            Otherwise, we will start with initial guess, then randomly pick 4 sets of HPs and initialize a probobalistic search
        """
        if init_vals=='None':
            trials = Trials()
            random=8
        else:
            trials = generate_trials_to_calculate(init_vals)
            random=4
        """Details for minimization """    
        best = fmin(fn=objective,
                    space=space,
                    algo=partial(tpe.suggest,n_startup_jobs=random),
                    max_evals=max_evals,
                    timeout=timeout, 
                    loss_threshold=loss_threshold,
                    trials=trials, 
                    show_progressbar=False)
        

        best_space=space_eval(space, best)
        

     
        RF = RandomForestRegressor(**{'max_depth': int(best['max_depth']),
                                      'n_estimators': int(best['n_estimators']),
                                      'random_state': int(best['random_state']),
                                      'min_samples_split': int(best['min_samples_split'])} ).fit(self.X_train, self.Y_train) 
    
        
        y_fit = RF.predict(self.X_train)
        y_predict = RF.predict(self.X_test)

        text=str("Time taken for RF: " + str(datetime.now() - startTime) ) 
        if verbose==1:
            text=str("Time taken for RF: " + str(datetime.now() - startTime) ) + str(modeling.verbose_output(self,  RF,  y_fit, y_predict))
  
        return  y_fit, y_predict, trials.trials, best_space, RF, text
 


  
    
    ####### later   def XGB(self,...):
        


               


