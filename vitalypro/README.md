**Prediction of Electronic Band-structure in Ultra-thin SiGe Superlattices with Tree-Based Machine Learning Models**

**Overview**

The idea of this task is to train an ML model to predict electronic band-structure of a SiGe superlattice with any kind of interfacial disorder, external strain, and composition*.   
(**please see the list of assumptions that will be used in this work*)

Superlattices are periodic structures that contain layers of different materials. The interfaces significantly impact the electronic and phonon transport and make superlattices ideal candidates as thermoelectric materials. The efficiency or figure of merit of thermoelectric material is calculated as:

*zT=(S^2/ρk)T*,

where *zT* is the dimensionless figure of merit, S, ρ, k  are the Seebeck coefficient, electrical resistivity, and thermal conductivity. In order to improve the efficiency of thermoelectric material, we need to increase the Seebeck coefficient and electrical conductivity and at the same time decrease the thermal conductivity.  The former two parameters can be easily and quickly computed from the energy bands using Boltzmann equations (see https://arxiv.org/pdf/cond-mat/0602203.pdf eq. 12-16). The most time-consuming part in the computation of S is the bandstructure calculation that can be significantly speed up with ML modeling. 
 

**Assumptions**

This task is very interesting, but it’s also complicated, and it would take a significant amount of time to complete. For this reason, I will make some assumptions that would simplify this problem:
1. I will only study superlattices that have a constant number of atoms to avoid band structure folding. As example, the electronic band structure unfolding can be done with GPAW (https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/electronic/unfold/unfold.html)
2. I will use a smaller number of atoms in the cell to speed up calculations (8 atom cells)
3. For the transport we only care about the bands that are close to the Fermi level. Therefore, I will only predict 2 valence and 2 conduction bands. 
4. In this work, I will only consider superlattices with the ideal interfaces. It would be interesting to investigate disordered interfaces and introduce some defects far away from the interface as well in the future work. This would require studying larger cells and introduce Voronoi tessellations to correctly describe the neighboring atoms. In this work, however, I will study the effect of composition and external strain only. 
5. I will predict bandstructure along a short path (e.g. Gamma -Z) with a small number of K-points
6. I will use the PBE exchange-correlation functional for the DFT calculations. The band gap will be underestimated but the bands shape should not be affected significantly. 
7. To build a good ML model we need to provide both global and local properties of the system. Global properties are lattice constants, compositions, number of atoms (if not a constant) etc. Local properties are such as local strain, nearest neighbor etc. In this work I will only use global properties to build an ML model  


**Project details**

1. **Step 1: DFT calculations** 

In the first step I performed DFT calculations of ultra-thin SixGe1-x superlattices at different external strain (x=1,2,3 is the number of monolayers). I will provide more details on this step during our Monday call.   

2. **Step 2: Develop a module for ML modeling**

This module consist of 3 classes: *read_data*, *data_preparation*, and *modeling*. The first class is written to get all necessary data for modeling from QE output files. This includes lattice constants, compositions, fermi level (read from xlm file) and band energies (read from bands.x output). *data_preparation* module can be used to split the data into training / testing sets. Options with manual or random selection are included. The final module consist of two tree-based ML models: Random forest and XGboost. Random forest suppots multioutput regression, however, for XGboost I had to use sklearn MultiOutputRegressor to make predictions (see https://github.com/dmlc/xgboost/issues/2087). The modeling is perfromed with automatic hyperparameter search using Bayesian optimization. To reduce possibility of overfitting, an option for n-fold cross validation is added. More details will be provided on Monday 
 
3. **Step 3: Predict a SL bandstructure with the ML model**

Finaly, when the structures are calculated, the data is extracted and prepared, we test the model performance. In the figure below, I present the final results for two SiGe SL bandstructures predicted with tree-based ML models. The black curves show the bands calculated with the DFT method, the predictions with random forest and XGboost models are shown with the red and blue dots, respectively. In general, both methods provide with good predictions, however, random forest model outperfroms the XGBoost one. I believe this mostly happens due to inability for XGBoost to deal with multioutputs. We will discuss this on Monday as well.  

![ml model performance](https://user-images.githubusercontent.com/64281595/139505384-c190b36a-62a3-48ee-9360-b1d56e0efeea.png)


