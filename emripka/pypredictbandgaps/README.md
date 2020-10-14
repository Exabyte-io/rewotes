# pyPredictBandGaps

* A package to predict a material's bandgap for a given set of stoichiometric and crystallographic properties. 
    * Relies on the Materials Project query funcionality for training data. 

## Dependencies
* **pymatgen**
    * follow installation instructions as suggested [here](https://pymatgen.org/index.html). 
    * Add your Materials Project API key: ``pmg config --add PMG_MAPI_KEY <USER_API_KEY>``

## Installation
1. Clone this repository: ``git clone https://github.com/emripka/rewotes.git``
2. Build the package inside the ``pypredictbandgaps`` outer directory: ``python3 setup.py build``
3. Install the package inside the same directory: ``python3 setup.py install``

## Usage 

1. Create a Material class for each material of interest 
    * Each object houses the novel material's stoichiometric and crystallographic properties
2. Use the BandGapPredicitons class to predict the band gap of each material. 
    * Feed in a list of Material objects, as well as a choice in training model. 

```python
# input materials
Si2Ge3 = Material(formula="Si2Ge3",spacegroup="F-43m",a=4.12,b=4.12,c=4.12,alpha=59.99,beta=59.99,gamma=59.99)
Si82Ge17 = Material(formula="Si82Ge17",spacegroup="F-43m",a=3.98,b=3.98,c=3.98,alpha=59.99,beta=59.99,gamma=59.99)

# train the model and predict the bandgaps
materials_list = [Cd2Se3, Si82Ge17]
band_gap_predictions = BandGapPredictions(materials_list, model_type="random_forest")
```
