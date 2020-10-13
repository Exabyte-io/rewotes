# pyPredictBandGaps

A package to predict a material's bandgap for a given set of stoichiometric and crystallographic properties. 

## Dependencies
**pymatgen**
    - follow installation instructions as suggested [here](https://pymatgen.org/index.html). 
    - Add your Materials Project API: ``pmg config --add PMG_MAPI_KEY <USER_API_KEY>``

## Installation
1. Clone this repository: ``git clone https://github.com/emripka/rewotes.git``
2. Build the package inside the ``pypredictbandgaps`` outer directory: ``python3 setup.py build``
3. Install the package inside the same directory: ``python3 setup.py install``

## Usage 

There are two distinct ways to use the package:

1. User inputs a novel material's stoichiometric and crystallographic properties, and the database of information at the Materials Project is used to train the model on the chemical system to which the material belongs. 

```python
# input materials
SiGe = Material(formula="SiGe",a=3.95,b=3.95,c=3.95,alpha=59.99,beta=59.99,gamma=59.99,volume=43.7409)
Si2Ge3 = Material(formula="Si2Ge3",a=4.25,b=4.25,c=4.25,alpha=59.99,beta=59.99,gamma=59.99,volume=45.222)

# train the model and predict the bandgaps
materials = [SiGe, Si2Ge3]
band_gap_predictions = BandGapPredictions(materials)
```

2. Alternativel,y the user inputs a novel material's stoichiometric and crystallographic properties, along with a input training data.  

```python
# input materials
SiGe = Material(formula="SiGe",a=3.95,b=3.95,c=3.95,alpha=59.99,beta=59.99,gamma=59.99,volume=43.7409)
Si2Ge3 = Material(formula="Si2Ge3",a=4.25,b=4.25,c=4.25,alpha=59.99,beta=59.99,gamma=59.99,volume=45.222)

# create training data
Si_training_data = TrainingData("Si",band_gap=0.514,a=3.95,b=3.95,c=3.95,alpha=59.99,beta=59.99,gamma=59.99,volume=43.7409)
Ge_training_data = TrainingData("Ge",band_gap=0.67,a=4.25,b=4.25,c=4.25,alpha=59.99,beta=59.99,gamma=59.99,volume=45.222)

# train the model and predict the bandgaps using the input training data
materials = [SiGe, Si2Ge3]
training_data = [Si_training_data, Ge_training_data]
band_gap_predictions = BandGapPredictions(materials, training_data)
```
