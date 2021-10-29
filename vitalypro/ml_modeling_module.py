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
1.  class for data (energy bands, cell geometry etc) extraction from QE output files
2.  class for data preparation for ML modeling
3.  class for ML modeling with RF and XGB 

""" 
class read_data:
    
    """ 
    In this module we read data from QE output and build a dataset for ML modeling 
    1. read energy band data from output produced by bands.x 
    2. read structure parameters and fermi level from xlm file
    In future work this module can be improved to read data from VASP and other modeling software 
    """
    def __init__(self, file):
        self.file = file
        


    def read_bands_ml(self, mode='Fermi', fermi=0, bands_num=2, npl=10):
        """
        This function reads energy bands from QE output file produced by band.x. 
        I found this code here https://github.com/yyyu200/QEbandplot/blob/master/pw_band_plot.py
        Modifications I made to the code:
            - bands get adjusted with respect to the fermi level
            - print only first n valence and n conduction bands above/below the fermi level 
        :mode - 'Fermi': user provides fermi level and bands get adjusted with respect to this level     
              - 'Auto' : the fermi level will be automatically placed at the top  of the conduction zone. 
              this option can be useful is the fermi level is unknown or calculated on a large mesh (my case)
        :fermi - fermi level (eV) (will be ignored if mode='Auto' is used.  
        :bands_num - number of valence/conduction bands near fermi level to print 
                
        """
        # feig : filband in bands.x input file
        # npl : number per line, 10 for bands.x, 6 for phonon
        import re
        f=open(self.file,'r')
        lines = f.readlines()
    
        header = lines[0].strip()
        line = header.strip('\n')
        shape = re.split('[,=/]', line)
        nbnd = int(shape[1])
        nks = int(shape[3])
        eig = np.zeros((nks, nbnd+1), dtype=np.float32)
    
        dividend = nbnd
        divisor = npl
        div = nbnd // npl + 1 if nbnd % npl == 0 else nbnd // npl + 2 
        kinfo=[]
        for index, value in enumerate(lines[1:]):
            value = value.strip(' \n')
            quotient = index // div
            remainder = index % div
    
            if remainder == 0:
                kinfo.append(value)
            else:
                value = re.split('[ ]+', value)
                a = (remainder - 1) * npl
                b = a + len(value)
                eig[quotient][a:b] = value
    
        f.close()
        if mode=='Fermi':
            eig = eig -fermi
            cond_bands=[]
            val_bands=[]
            for i in range(eig.shape[0]):
                cond_bands=np.append(cond_bands, eig[i][eig[i]>0][:bands_num])    
                val_bands =np.append(val_bands,eig[i][eig[i]<0][-bands_num-1:-1])
        
            cond_bands=cond_bands.reshape(11,bands_num).swapaxes(0,1)
            val_bands=val_bands.reshape(11,bands_num).swapaxes(0,1)
            bands=np.append(val_bands,cond_bands)
        elif mode=='Auto':
            """
            For a general case, i would need to read QE output file or pseudopotential file to find out
            how many valence electrons each atoms has and how many atoms of each type I have to count the number 
            of conduction bands. It can be done in the future work. In this work I only study SiGe systems. Both Si and Ge have 
            4 val electrons. All my structures contain 8 atoms, so 4*8/2=16 conduction bands
            number_cnb = 16
            """
            number_cnb = 16
            cond_bands=[]
            val_bands=[]
            for i in range(eig.shape[0]):
                cond_bands=np.append(cond_bands, eig[i][number_cnb-bands_num:number_cnb])    
                val_bands =np.append(val_bands,eig[i][number_cnb: number_cnb+bands_num])
            fermi= np.max(cond_bands)
            cond_bands=cond_bands-fermi
            val_bands=val_bands-fermi
            cond_bands=cond_bands.reshape(11,bands_num).swapaxes(0,1)
            val_bands=val_bands.reshape(11,bands_num).swapaxes(0,1)
            bands=np.append(val_bands,cond_bands)
            
        return bands

    def read_lattice_parameters(self):
        """read lattice constants from xlm file"""
        import xml.etree.ElementTree as ET
        from collections import Counter
        tree = ET.parse(self.file)
        root = tree.getroot()
        a = root[2][2][1][0].text
        c = root[2][2][1][2].text
        a = [float(value) for value in a.split(' ')]
        c = [float(value) for value in c.split(' ')]
        
        #read  composition SixGe1-x
        atom=[]
        for i in range(8):
            atom=np.append(atom,root[2][2][0][i].attrib['name'])    
        composition = int(Counter(atom)['Si'])/8
        
        return a[0], c[2], composition
    
    def read_fermi_levle(self):
        import xml.etree.ElementTree as ET
        tree = ET.parse(self.file)
        root = tree.getroot()
        """
        read fermi level from xlm file
        convert from Ry to eV and multiply by 2!!! 
        I am using qe 6.4 which has a bug for fermi level:
    
        Fermi energy incorrectly written to xml file in 'bands' calculation
        (did not affect results, just Fermi energy position in band plotting)
        see for more info https://gitlab.com/QEF/q-e/-/releases
        
        """
        fermi_level=round(float(root[3][9][7].text)*2*13.605,2)   
        return fermi_level

class data_preparation:
    
    """ 
    In this module we prepare data collected from QE output files for modeling. 
    The following functions split the data into training / testing sets.
    """
    def __init__(self, X, Y):
        self.X = X
        self.Y = Y
           
    def random_split(self,  train_split=0.8, seed=42): 
         """ 
         random_split function allows to randomly select training and testing sets from the data. 
         :train_split - shows the ratio of train set to to the whole set in the split. Can be from 0 to 1. Defauls is 0.8
         :seed is a random state for reproducibility  
         """
         from sklearn.model_selection import train_test_split
         X_train, X_test, Y_train, Y_test = train_test_split(self.X, self.Y, test_size=float(1-train_split), random_state=seed)
         return X_train, X_test, Y_train, Y_test
    
    def select_train_test(self, test_row_number=[5,6]):
        """ 
        select_train_test function allows to manually specify training and testing sets from the data. 
        :test_row_number is the number of rows that will be used for testing. The rest is training data
    
        """
        X_test = self.X[test_row_number]
        Y_test = self.Y[test_row_number]
        X_train = np.delete(self.X, test_row_number, 0)
        Y_train = np.delete(self.Y, test_row_number, 0)
        return X_train, X_test, Y_train, Y_test
    
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
        
        :max_evals - number of max evaluation in hyperparameter search algorithm. The search will stop when the number of max itteration is reached. Can be any integer number. Default value is 40, 
        :timeout - search time. The search will stop when the time is out. Can be any integer number. Default value is 600 sec  
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
           More hyperparameters to build the space can be added in the future work """    
        
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
 


  
    def XGB(self,  max_evals=40, timeout=600,  loss_threshold=0.5, init_vals='None' , n_jobs=10, cv=5, verbose=1):
        

        """
        XGBoost does not support multi-output predictions. In this example I wanted to demonstrate one of the possible solutions 
        for energy bands prediction with ML methods like XGB. 
        For XGB model we wil use the same strategy as we used for Random forest model. In addtion we will apply sklearn  MultiOutputRegressor 
        for multioutput predictions. This method is VERY SLOW...
        
        """
        

        import xgboost as xgb
        from sklearn.multioutput import MultiOutputRegressor

        startTime = datetime.now()
        
        def objective(space):
    
            clf = MultiOutputRegressor(xgb.XGBRegressor(n_estimators =int(space['n_estimators']),                           
                                   learning_rate = .25,
                                   max_depth = int(space['max_depth']),
                                   subsample = space['x_subsample'],
                                   alpha = space['alpha'],
                                   reg_lambda = space['x_reg_lambda'],
                                   min_split_loss = space['min_split_loss'],
                                   objective='reg:squarederror'))
    
    
            mse_scr = cross_val_score(clf, self.X_train, self.Y_train, scoring='neg_root_mean_squared_error', cv=5,verbose=0,n_jobs=n_jobs).mean()

            return {'loss': -mse_scr, 'status': STATUS_OK }
    
    
        space ={
                'max_depth':  hp.quniform('max_depth', 1, 18, 1),
                'x_subsample': hp.quniform ('x_subsample', 0.2, 1, 0.02),
                'alpha' : hp.quniform ('alpha', 0.0, 0.6, 0.02),
                'x_reg_lambda' : hp.quniform ('x_reg_lambda', 0, 1, 0.02),
                 'min_split_loss' : hp.quniform ('min_split_loss', 0, 8, 0.02),
                'n_estimators' : hp.quniform('n_estimators', 40, 800, 10)
            }
        
        if init_vals=='None':
            trials = Trials()
            random=8
        else:
            trials = generate_trials_to_calculate(init_vals)
            random=4
            

        best = fmin(fn=objective,
                    space=space,
                    algo=partial(tpe.suggest,n_startup_jobs=random),
                    max_evals=max_evals,
                    timeout=timeout, 
                    loss_threshold=loss_threshold,
                    trials=trials, 
                    show_progressbar=False)

        best_space=space_eval(space, best)

    
    
        regressorXGB = MultiOutputRegressor(xgb.XGBRegressor(**{'max_depth': int(best_space['max_depth']),
                                           'n_estimators': int(best_space['n_estimators']),
                                           'x_subsample': float(best_space['x_subsample']),
                                           'alpha': float(best_space['alpha']),
                                           'x_reg_lambda': float(best_space['x_reg_lambda']),
                                           'min_split_loss': float(best_space['min_split_loss'])  },                                     
                                          objective='reg:squarederror')).fit(self.X_train, self.Y_train) 
                                           
        y_fit = regressorXGB.predict(self.X_train)
        y_predict = regressorXGB.predict(self.X_test)

        text=str("Time taken for XGB: " + str(datetime.now() - startTime) ) 
        if verbose==1:
            text=str("Time taken for XGB: " + str(datetime.now() - startTime) ) + str(modeling.verbose_output(self,  regressorXGB,  y_fit, y_predict))
  
        return  y_fit, y_predict, trials.trials, best_space, regressorXGB, text
        


               


