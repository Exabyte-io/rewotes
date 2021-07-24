# Discussion

We have provided code to query structures and relevant properties for the Materials Project API. In order to reduce the
number of requests we save the structures using a utility from matminer which ensures that pymatgen objects are correctly
serialised.

The first set of modelling code we have produced - `composition.py` -- can be used to train compositon based machine learning models.
We use the magpie feature set [1] which we construct using `matminer` and a `RandomForestRegressor` from `sklearn`. We augment the `sklearn`
implementation to provide a crude uncertainty estimate from the leaf impurity in the trees, this code is stored in `utils.py` as it is
reused in other scripts. Generally speaking this code could be improved by carrying out hyperparameter searching. However, the added value
from this in this case is low due to the fact that the composition-based model is ill-suited to the task presented of investigating
the bandgap of the Si-Ge chemical system. This can be seen by looking at the two figures "SiGe-composition-combined.pdf" and
"SiGe-composition-non-metals.pdf" (produced using `explore-composition.py`) where we see that the model trained on all the data performs very poorly. Whilst the model trained
on the non-metal data only seemingly offers better performance this improvement is not reliable as evidenced by the large uncertainty
in the model predictions.

Due to the limitations of a composition-based approach we next look at structure-based models. Various models have been reported in the
literature for predicting the bandgap of materials given their structures. Structure-based machine learning models typically operate
on local atomic environments performing pooling to make a single prediction for an the whole material. We make use of the `MEGNet` model
due to it's tight integration with `pymatgen` and other tools from the Materials Project team. Generally `MEGNet` is easy to use and provides
several pretrained sets of model weights. One limitation of `MEGNet` is that the code for generating the neighbour lists and distances
(i.e. the graphs on which the graph neural network operates) is handled inside `pymatgen`. This means that the `MEGNet` codebase is not well
suited to the prediction of forces as the co-ordinates are not present in the computation graph meaning that the forces cannot be computed
by automatic differentation of the energy with respect to the input co-ordinates. Other codebases do implement this functionality but typically
the models are implemented in a highly opinionate manner with a view to application in ML-FF/MD applications.As such they are not appropraite
for the band gap prediction task given.

Due to computational resources I have opted to use the pretrained band gap model provided by `MEGNet` for futher exploration of how structure
based ML models perform in the Si-Ge chemical system rather than training a model on all the data from scratch. The Si-Ge chemcial system is
as both Si and Ge adopt the same ground state structure and there are no other stable structures in the phase diagram at the DFT level using a
PBE functional. The experimental band gap (as measured and parameterised in [2]) shows a cross over in behaviour at Si0.15Ge0.85. We want to
test whether the limited data in the Si-Ge system available in the Materials Project is sufficient to recover this beahviour. The expectation
is that the structure-based model will fail to reflect the experimental reality as it is trained on a sparse set of DFT calculated band gaps
that do not reflect experiment themselves. Most notably the ground state for Ge in the Materials Project is calculated to be metallic.

As structure-based neural networks are extensible and pool over local environments we should in principle be able to estimate a range of possible
bandgaps at different compositions by predicting the band gap of different random ordered-disordered superstructures. We could then make an estimate of
bandgap by taking a boltzman weighted mean over the different distinct superstructures -- we would likewise estimate these energies using a pretrained
`MEGNet` model. I attempted to implement this investigation by using the `EnumerateStructureTransformation` implemented in `pymatgen` however I encountered
errors whilst compiling the `enumlib` package used by the transformation class. I was unable as of yet to resolve this problem in the time spent
on the rewote task -- the most likely path to a solution would be to contact the `enumlib` developers.

In general the task of predicting the band gaps of materials using a machine learning model is a multi-faceted problem. In this exercise I have attempted
to highlight some aspects of this. It is possible to build predictive models for the band gap that achieve high coefficients of determination
for both composition and structure based approaches. Such models are sufficiently accurate (MAE ~0.3-0.5 eV) for applications in materials discovery.
However, for industrial relevant tasks looking as specific systems the general purpose models trained on chemically diverse data sets are unlikely
to be useful as this level of error is beyond industrially acceptable tolerances. In contrast to other problems I view it as a more bespoke task that
will be challenging to automate due to its open ended nature.

# Time breakdown

* 2 hours were spent reviewing literature material for Si-Ge and choosing which Machine Learning models to use.
* 2 hours were spent exploring test models and ideas in a google colab notebook before the indermediate meeting.
* 2 hours were spent writing/transfering these ideas into the (working) scripts provided in this repo.
* 1 and a half hours were spent trying to debug the `enumlib` issue that is blocking the planned deliverables with respect to the structure-based exploration

Total: 7 hrs 30 mins

# References

[1] https://www.nature.com/articles/npjcompumats201628
[2] https://journals.aps.org/prb/abstract/10.1103/PhysRevB.40.5683
