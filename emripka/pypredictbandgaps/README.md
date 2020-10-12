# pyPredictBandGaps

A package to predict a material's bandgap for a given set of stoichiometric and crystallographic properties.

## Installation
1. Clone this repository: ``git clone https://github.com/emripka/rewotes.git``
2. Build the package inside the ``pypredictbandgaps`` outer directory: ``python3 setup.py build``
3. Install the package inside the same directory: ``python3 setup.py install``

## Usage 

There are two distinct ways to use the package:

1. User inputs a novel material's stoichiometric and crystallographic properties, and the built-in database of training data is used to train the model. 

```python
# input materials
SiGe = Material(formula="SiGe",density=1.234,crystal_structure="cubic")
Si2Ge3 = Material(formula="Si2Ge3",density=2.135,crystal_structure="cubic")

# train the model and predict the bandgaps
materials_list = [SiGe, Si2Ge3]
band_gap_predictions = BandGapPredictions(materials_list)
```

2. User inputs a novel material's stoichiometric and crystallographic properties, as well as training data, and choose to use only this new training data to train the model, or use it in conjunction with the database. 

```python
# input materials
SiGe = Material(formula="SiGe",density=1.234,crystal_structure="cubic")
Si2Ge3 = Material(formula="Si2Ge3",density=2.135,crystal_structure="cubic")

# create training data
Si_training_data = TrainingData("Si",band_gap=0.514,density=2.282,crystal_structure="hexagonal")
Ge_training_data = TrainingData("Ge",band_gap=0.67,density=4.445,crystal_structure="cubic")

# train the model and predict the bandgaps
# using only the input training data
materials_list = [SiGe, Si2Ge3]
input_training_data = [Si_training_data, Ge_Trainig_data]
band_gap_predictions = BandGapPredictions(materials_list,input_training_data,use_database_data=False)
```
