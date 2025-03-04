# ReWoTes: ML Property Predict

Alec Glisman

## Overview

This directory contains files for the ML Property Predict project for Mat3ra.com.

Input data is accessed from the Materials Project and the data is cleaned into Pandas Dataframes inside `data/data_load.py`.
I chose to download all materials with a bandgap of less than 10 eV from the Materials Project and parsed all data related to the crystallographic and stoichiometric properties.
Categorical data is converted to numeric data using one-hot encoding and the data is then scaled using `sklearn.preprocessing.StandardScaler`.
The input data source to the machine learning model can be augmented with additional Materials Project data with the `MaterialData` init method and external data can also be merged using its respective `add_data_columns` method.
The cleaned data is archived using Pandas in conjunction with HDF5 to lower runtime costs for model development.

I chose to pursue two machine-learning architectures: XGBoost and feed-forward, fully connected, neural networks.
XGBoost generally performs better than neural networks when the data set is not large, and XGBoost is also much faster to train.
Neural networks were included for their superior expressivity and serve as a useful comparison to XGBoost.
In both cases, I employed `KFold` and `RandomizedSearchCV` from `scikit-learn` to cross-validate and select hyperparameters, respectively.

The best XGBoost Regressor that I trained is saved during runtime under the `models` directory and has a testing sample MSE of 0.646 eV.
Similarly, the best fully connected neural network I trained is saved during runtime under the `models` directory and has a testing sample MSE of 0.817 eV.
The seed used is provided in `main.py` for reproducibility.

Areas for future work include:

1. Stratified sampling for test/train split or cross-validation to make sure different space groups are represented properly in each subset.
2. Explore the use of feed-forward neural networks and experiment with architecture, drop-out, and regularization to optimize the performance. Additionally, increase the epochs from 40. I used 40 due to computational constraints, but the loss was still noticeably shrinking.
3. Addition of more data from the Materials Project to lower the inductive bias of the models.
4. Attempt transfer-learning of these models and fine-tune to more specific databases, such as silicon semiconductors.

## Usage

A Conda environment file has been provided (`requirements.yml`) to set up a Python environment called `ml-band-gaps` with the following command

```[bash]
$ conda env create -f requirements.yml
```

The overall project can then be run with

```[bash]
$ python main.py
```

Unit tests can be run with pytest as

```[bash]
$ pytest tests
```

Data ingested is cached to the `data` directory, and machine-learning models are cached to the `models` directory.
Each of these directories is created automatically as part of the main script.

Note that the data is sourced from the Materials Project, which requires an API key to access it.
I have added my API key to the `.gitignore` for security reasons, so users will need to generate their own and add it to an `api_key.txt` file.

## Requirements

- suggest the bandgap values for a set of materials designated by their crystallographic and stoichiometric properties
- the code shall be written in a way that can facilitate easy addition of other characteristics extracted from simulations (forces, pressures, phonon frequencies etc.)

## Expectations

- the code shall be able to suggest realistic values for slightly modified geometry sets - e.g. trained on Si and Ge it should suggest the value of bandgap for Si49Ge51 to be between those of Si and Ge
- modular and object-oriented implementation
- commit early and often - at least once per 24 hours
