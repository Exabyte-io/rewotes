# ML-Band-Gaps

## Data

* The data source we are using for this work is the Materails Project. MP uses GGA (PBE) planewave DFT implemented in VASP. In total they have conducted energy calculations for ~130k Materials.
* Machine Learning models are only able to approximate the quality of the data used to train them. If we use DFT data to train we can only approximate the accuracy of the level of theory of the training data calculations. DFT band gaps predictions are notorious for underestimating the True band gap particularly GGA functions, higher levels of theory i.e. mGGA, Hybrid functional, GW produce better quality predictions at significantly increased cost.
* We want to avoid covariate shift in our ML model i.e. changes in DFT settings not reflected in model inputs. MP reports band gaps in some instances when they haven't run their "high quality" band structure workflow -- Only use data satisfying where this workflow has been run.
* In the Si-Ge phase diagram there are no stable phases apart from the pure phases in MP. DFT also gives the ground state structure of Ge as a metal.

## Frameworks

* As we're using MP as a data source we will use frameworks built around pymatgen. This means using Matminer for shallow learning applications and MEGNet for Structure-based deep learning applications.
* We will first train a composition based model on all the data queried and show it's predictions in the Si-Ge chemical system. For this we will use the Magpie descriptor set.
* MEGNet provides a pre-trained bandstructure model from 2018. For the sake of this task we will use this model for downstream investigations but provide code that would enable such a model to be trained from scratch.

## Deliverables

- [x] Code to query relevant data from the Materials Project.
- [x] Code to train a Random Forest Regression model (sklearn) using the Magpie descriptor set (Matminer).
- [ ] Code to train a MEGNet regression model (tensorflow, megnet).
- [x] Code to generate predictions over the entire compositon space for compositon based models.
- [ ] Code to generate supercells in the Si-Ge phase space to explore the variations in predictions of structure-based models.
- [ ] Dicussion document.
