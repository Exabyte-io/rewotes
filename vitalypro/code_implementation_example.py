#!/usr/bin/env python
# coding: utf-8

# In[1]:


#load some modules  
from ml_modeling_module import *
from os import walk
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# # first let's plot some sige bandstructures calculated with DFT

# In[2]:


### EXAMPLE OF READING DATA WITH MODE="FERMI" ###

"""
read QE output file with bands for 0% and 5% strain 
reshape data: 8 bands 11 k-points
fermi level is manually specified 
"""
data=read_data('bands/si2ge2_0.dat')
si2ge2_0=data.read_bands_ml(mode='Fermi',fermi=6.18, bands_num=4)
si2ge2_0=si2ge2_0.reshape(8,11)

data=read_data('bands/si2ge2_5.dat')
si2ge2_5=data.read_bands_ml(mode='Fermi',fermi=5.60, bands_num=4)
si2ge2_5=si2ge2_5.reshape(8,11)

fig = make_subplots(rows=1, cols=2,subplot_titles=("Si2Ge strain 0%", "Si2Ge2 strain 5%"))

for i in range(si2ge2_0.shape[0]):
    fig.add_trace(go.Scatter(x=np.arange(1,12,1), y=si2ge2_0[i],
                        mode='lines', marker=dict(color='black'), showlegend=False, 
                        name='markers'), col=1, row=1)  
for i in range(si2ge2_5.shape[0]):
    fig.add_trace(go.Scatter(x=np.arange(1,12,1), y=si2ge2_5[i],
                        mode='lines', marker=dict(color='black'), showlegend=False, 
                        name='markers'), col=2, row=1) 

    
fig.add_trace(go.Scatter(x=[1,11], y=[0,0],
                        mode='lines', line=dict(color= "RoyalBlue",width=3, dash='dot'), showlegend=False, 
                        name='markers'), col=1, row=1) 
fig.add_trace(go.Scatter(x=[1,11], y=[0,0],
                        mode='lines', line=dict(color= "RoyalBlue",width=3, dash='dot'), showlegend=False, 
                        name='markers'), col=2, row=1) 
fig.update_xaxes(showgrid=True,  
                 gridwidth=0.5, 
                 gridcolor='lightgrey', 
                 ticks='inside' ,                        
                 linecolor= 'black',                      
                 linewidth= 1,
                 mirror= True,
                 ticktext=['G', 'Z'],
                 tickvals=[1,11],
                 range=[1,11],
                 )
fig.update_yaxes(title='Energy, eV', 
                 showgrid=True,  
                 gridwidth=0.5, 
                 gridcolor='lightgrey', 
                 ticks='inside' ,                        
                 linecolor= 'black',
                 linewidth= 1,
                 mirror= True,
                )

fig.update_layout(width=900, height=650, plot_bgcolor='white')
fig.show()

"""
I used only 11 k-points to produce these bandstructures,
therefore, the bandstructures look a bit ugly
We can see how conduction bands shift with applied strain
"""


# In[3]:


### EXAMPLE OF READING DAT WITH MODE="Auto" ###

"""
read QE output file with bands for 0% and 5% strain and different composition
reshape data: 8 bands 11 k-points
In this case the fermi level is not specified,
it is automatically placed at the top of the valence zone
"""

data=read_data('bands/si1ge3_0.dat')
si1ge3_0=data.read_bands_ml(mode='Auto', bands_num=4)
si1ge3_0=si1ge3_0.reshape(8,11)

data=read_data('bands/si1ge3_5.dat')
si1ge3_5=data.read_bands_ml(mode='Auto', bands_num=4)
si1ge3_5=si1ge3_5.reshape(8,11)

data=read_data('bands/si3ge1_0.dat')
si3ge1_0=data.read_bands_ml(mode='Auto', bands_num=4)
si3ge1_0=si3ge1_0.reshape(8,11)

data=read_data('bands/si3ge1_5.dat')
si3ge1_5=data.read_bands_ml(mode='Auto', bands_num=4)
si3ge1_5=si3ge1_5.reshape(8,11)


fig = make_subplots(rows=2, cols=2,
                    subplot_titles=("Si1Ge3 strain 0%", " Si1G3 strain 5%", 'S3Ge1 strain 0%', 'Si3Ge1 strain 5%'),
                   vertical_spacing=0.04)

for i in range(si3ge1_0.shape[0]):
    fig.add_trace(go.Scatter(x=np.arange(1,12,1), y=si3ge1_0[i],
                        mode='lines', marker=dict(color='black'), showlegend=False, 
                        name='markers'), col=1, row=1)  
for i in range(si3ge1_5.shape[0]):
    fig.add_trace(go.Scatter(x=np.arange(1,12,1), y=si3ge1_5[i],
                        mode='lines', marker=dict(color='black'), showlegend=False, 
                        name='markers'), col=2, row=1) 
for i in range(si1ge3_0.shape[0]):
    fig.add_trace(go.Scatter(x=np.arange(1,12,1), y=si1ge3_0[i],
                        mode='lines', marker=dict(color='black'), showlegend=False, 
                        name='markers'), col=1, row=2) 
for i in range(si1ge3_5.shape[0]):
    fig.add_trace(go.Scatter(x=np.arange(1,12,1), y=si1ge3_5[i],
                        mode='lines', marker=dict(color='black'), showlegend=False, 
                        name='markers'), col=2, row=2) 
    
fig.add_trace(go.Scatter(x=[1,11], y=[0,0],
                        mode='lines', line=dict(color= "RoyalBlue",width=3, dash='dot'), showlegend=False, 
                        name='markers'), col=1, row=1) 
fig.add_trace(go.Scatter(x=[1,11], y=[0,0],
                        mode='lines', line=dict(color= "RoyalBlue",width=3, dash='dot'), showlegend=False, 
                        name='markers'), col=2, row=1) 
fig.add_trace(go.Scatter(x=[1,11], y=[0,0],
                        mode='lines', line=dict(color= "RoyalBlue",width=3, dash='dot'), showlegend=False, 
                        name='markers'), col=1, row=2) 
fig.add_trace(go.Scatter(x=[1,11], y=[0,0],
                        mode='lines', line=dict(color= "RoyalBlue",width=3, dash='dot'), showlegend=False, 
                        name='markers'), col=2, row=2) 
fig.update_xaxes(showgrid=True,  
                 gridwidth=0.5, 
                 gridcolor='lightgrey', 
                 ticks='inside' ,                        
                 linecolor= 'black',                      
                 linewidth= 1,
                 mirror= True,
                 ticktext=['G', 'Z'],
                 tickvals=[1,11],
                 range=[1,11],
                 )
fig.update_yaxes(title='Energy, eV', 
                 showgrid=True,  
                 gridwidth=0.5, 
                 gridcolor='lightgrey', 
                 ticks='inside' ,                        
                 linecolor= 'black',
                 linewidth= 1,
                 mirror= True,
                )

fig.update_layout(width=900, height=1450, plot_bgcolor='white')
fig.show()

"""
I used only 11 k-points to produce these bandstructures,
therefore, the bandstructures look a bit ugly
We can see how conduction bands shift with applied strain
"""


# In[4]:


""" 
reading cell parameters from all QE output xlm files
"""

mypath='cell_param/'
filenames = next(walk(mypath), (None, None, []))[2]

lattice=[]
fermi=[]
for file in filenames:
    data=read_data(mypath+file)
    (a, c, composition) = data.read_lattice_parameters()
    lattice.append ([a, c, composition])
# Fermi energy could be obtained with:    
#    fermi_level = data.read_fermi_levle()
#    fermi.append(fermi_level)
    
X=lattice=np.array(lattice)

X


# In[36]:


filenames


# In[5]:


''' 
reading bands from all QE output files 
(output files are produced with band.x module)
'''

mypath='bands/'
filenames = next(walk(mypath), (None, None, []))[2]

y=[]
for i in range(len(filenames)):
    data=read_data(mypath+filenames[i])
    bands_ml=data.read_bands_ml(mode='Auto', bands_num=2)
    y= np.append(y, bands_ml)
y=y.reshape(len(filenames),44)
y.shape


# # prepare data for modeling: split for training and testig sets

# In[8]:


my_data=data_preparation(X,y)

#to manually specify testing set, use the following:
(X_train, X_test, Y_train, Y_test)=my_data.select_train_test(test_row_number=[12,16])
# for a random split use function:
#(X_train, X_test, Y_train, Y_test)=my_data.random_split(train_split=0.8, seed=42)

my_model=modeling(X_train, Y_train, X_test, Y_test, scaler='None')

#train Random forest model with No initial guess
(y_fit_RF, y_predict_RF, trials_RF, best_space_RF, regressorRF, text_RF)=my_model.RandomForest( max_evals=12, timeout=150*60,loss_threshold=0.01, cv=5, verbose=1)

print(text_RF)


#


# In[9]:


## MultiOutputRegressor is VERY SLOW!
## need more time to better understand MultiOutputRegressor

#train XGboost model with initial guess
best_space = {'alpha': 0.0, 'max_depth': 3.0, 'min_split_loss': 5.6, 'n_estimators': 100, 'x_reg_lambda': 1.00, 'x_subsample': 1.0}
(y_fit_XG, y_predict_XG, trials_XG, best_space_XG, regressorXGB, text_XG)=my_model.XGB( max_evals=12, timeout=150*60,loss_threshold=0.01,cv=5, verbose=1, init_vals=[best_space])
print(text_XG)


# # let's plot predicted bands vs calculated bands

# In[11]:


#reading calculated bands:
data=read_data('bands/si2ge2_relax.dat')
bands_si2ge2_relaxed=data.read_bands_ml(mode='Auto', bands_num=4)
bands_si2ge2_relaxed=bands_si2ge2_relaxed.reshape(8,11)

data=read_data('bands/si3ge1_3.dat')
bands_si3ge1_3=data.read_bands_ml(mode='Auto', bands_num=4)
bands_si3ge1_3=bands_si3ge1_3.reshape(8,11)

#predicted bands
bands_RF_1=y_predict_RF[0].reshape(4,11)
bands_RF_2=y_predict_RF[1].reshape(4,11)

bands_XG_1=y_predict_XG[0].reshape(4,11)
bands_XG_2=y_predict_XG[0].reshape(4,11)

fig = make_subplots(rows=1, cols=2,subplot_titles=("Si2Ge2 strain 2.5%", "Si3Ge1 strain 3%"))

for i in range(bands_si2ge2_relaxed.shape[0]):
    fig.add_trace(go.Scatter(x=np.arange(1,12,1), y=bands_si2ge2_relaxed[i],
                        mode='lines', marker=dict(color='black'), showlegend=False, 
                        name='markers'), row=1, col=1)
    
for i in range(bands_RF_1.shape[0]):
    fig.add_trace(go.Scatter(x=np.arange(1,12,1), y=bands_RF_1[i],
                        mode='markers', marker=dict(color='red'), showlegend=False, 
                        name='markers'), row=1, col=1)
for i in range(bands_XG_1.shape[0]):
    fig.add_trace(go.Scatter(x=np.arange(1,12,1), y=bands_XG_1[i],
                        mode='markers', marker=dict(color='blue'), showlegend=False, 
                        name='markers'), row=1, col=1)
for i in range(bands_si3ge1_3.shape[0]):
    fig.add_trace(go.Scatter(x=np.arange(1,12,1), y=bands_si3ge1_3[i],
                        mode='lines', marker=dict(color='black'), showlegend=False, 
                        name='markers'), row=1, col=2)
    
for i in range(bands_RF_2.shape[0]):
    fig.add_trace(go.Scatter(x=np.arange(1,12,1), y=bands_RF_2[i],
                        mode='markers', marker=dict(color='red'), showlegend=False, 
                        name='markers'), row=1, col=2)   
for i in range(bands_XG_2.shape[0]):
    fig.add_trace(go.Scatter(x=np.arange(1,12,1), y=bands_XG_2[i],
                        mode='markers', marker=dict(color='blue'), showlegend=False, 
                        name='markers'), row=1, col=2) 
    
    
fig.add_trace(go.Scatter(x=[1,11], y=[0,0],
                        mode='lines', line=dict(color= "RoyalBlue",width=3, dash='dot'), showlegend=False, 
                        name='markers'), col=1, row=1)    
fig.add_trace(go.Scatter(x=[1,11], y=[0,0],
                        mode='lines', line=dict(color= "RoyalBlue",width=3, dash='dot'), showlegend=False, 
                        name='markers'),col=2, row=1)  

fig.add_trace(go.Scatter(x=[12,13], y=[0,0],
                        mode='lines', line=dict(color= "black"), showlegend=True, name='DFT-PBE', 
                        ),col=1, row=1) 
fig.add_trace(go.Scatter(x=[12], y=[0],
                        mode='markers', line=dict(color= "red"), showlegend=True, name='RandomForest', 
                        ),col=1, row=1) 
fig.add_trace(go.Scatter(x=[12], y=[0],
                        mode='markers', line=dict(color= "blue"), showlegend=True, name='MultiOutputRegressor-XGBoost', 
                        ),col=1, row=1) 

fig.update_xaxes(showgrid=True,  
                 gridwidth=0.5, 
                 gridcolor='lightgrey', 
                 ticks='inside' ,                        
                 linecolor= 'black',                      
                 linewidth= 1,
                 mirror= True,
                 ticktext=['G', 'Z'],
                 tickvals=[1,11],
                 range=[1,11]
                 )
fig.update_yaxes(title='Energy, eV', 
                 showgrid=True,  
                 gridwidth=0.5, 
                 gridcolor='lightgrey', 
                 ticks='inside' ,                        
                 linecolor= 'black',
                 linewidth= 1,
                 mirror= True,
                )
fig.update_layout(title='Si2Ge2 relaxed', width=1000, height=650, plot_bgcolor='white')
fig.show()


# # example of solving reverse problem 

# In[25]:


# in this example we predict lattice constants for a given bandstructure (lattice constants are in bohr)

import GPyOpt
from GPyOpt.methods import BayesianOptimization
from sklearn.metrics import mean_squared_error

data=read_data('bands/si2ge2_relax.dat')
bands=data.read_bands_ml(mode='Auto', bands_num=2)


def objfunc(x):

    x1 = float(x[0,0])
    x2 = float(x[0,1])
    x3 = float(x[0,2])

    return  mean_squared_error(regressorRF.predict([[x1,x2,x3]]).reshape(1,44), bands.reshape(1,44))

maxiter = 40
domain =    [
             {'name': 'a_paral', 'type': 'continuous', 'domain': (10, 11)},
             {'name': 'a_perp',  'type': 'continuous', 'domain': (10, 11)},
             {'name': 'composition',  'type': 'discrete', 'domain': (0.25,0.5,0.75)}
            ]



Bopt = GPyOpt.methods.BayesianOptimization(objfunc, domain=domain)
Bopt.run_optimization(max_iter = maxiter)

print("Lattice parameters:"+str(Bopt.x_opt))    
print("MSE: "+str(Bopt.fx_opt)) 

"""
I will explain this part on Monaday 
"""


# # Example with random split 80/20

# In[13]:


# one more example to train model with a random split 80/20
my_data=data_preparation(X,y)
#random split could also be used
(X_train, X_test, Y_train, Y_test)=my_data.random_split(train_split=0.8, seed=42)

my_model=modeling(X_train, Y_train, X_test, Y_test, scaler='None')

#train Random forest model with No initial guess
(y_fit_RF, y_predict_RF, trials_RF, best_space_RF, regressorRF, text_RF)=my_model.RandomForest( max_evals=12, timeout=150*60,loss_threshold=0.01, cv=5, verbose=1)

print(text_RF)

# 20% of structures are predicted with a high accuracy

