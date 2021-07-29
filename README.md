# SiGe Properties Estimator (ML, Materials)

>Based on SiGe lattice parameters (flask server)

## Overview

The aim of this task is to create a python package that implements automatic prediction of various properties for a SiGe alloys, such as electronic band gaps, formation energy, bulk modulus, based on user's training data and test data. 

Target properties can vary and will be detrmined directly from a given training dataset by the estimator, i.e. per available data, avoiding not related features. 

### Usage

Minimum features required by the model are: lattice size in A, angles and lattice volume in A3, have to be present in the training dataset with following names:

- SiGe formula, as Si49Ge51 or Si, Ge alone
- a/A, b/A, c/A
- alpha, beta, gamma
- Volume/A3

Estimator will allocate rest of the features as target fingerprints for modeling, then will train models based on each of them one after another. ML model is based on SVR, suitable for smaller datsets (assumed as the most common for such modeling type). The trained models results are shown under route /model_train. Below are few suggested properties for the estimation:

- Band_Gap
- Bulk_modulus
- Thermal_expan_coeff
- Melting_point
- Thermal_conductivity
- Dielectric_const
- Optical_photon_energy
- Density
- Surface_microhardness
- Formation energy

Estimator have models tuning option under route /model_tune. Training datset size can be as low as 10 rows, up to 10,000 for reasonable wating times. Testing datset size can be as low as one row and can only include minimum required features. If no target fingerprints are provided in the test dataset, model will only predict the values, when the model score will be based on the split from the training dataset per target fingerprints found in the training dataset. 


