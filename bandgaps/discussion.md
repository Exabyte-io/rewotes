# Discussion

## Technical Discussion

Initial work was carried out in google colab notebook before being split out into separate scripts. In this instance I have not adopted a particularly object oriented setup due to the nature of the task.

* `query.py` - this script makes a query to the MaterialsAPI using the `MPRester` tool from `pymatgen`. The resulting data is placed into a `pandas` `DataFrame` and serialized using `io` tools from `matminer` that correctly serialize `MSONable` objects -- i.e. `pymatgen` objects such as `Structure` and `Composition`.

* `composition.py` - this script trains a `sklearn` `RandomForestRegressor` that takes the Magpie feature set [1] as input. The Magpie features are evaluated using `matminer`. We wrap the `RandomForestRegressor` implemented in `sklearn` to add in uncertainty estimation via leaf impurity - this has no impact on model training.

* `explore-composition.py` - this script plots the band gap as a function of composition in the Si-Ge chemical system.

* `structure.py` - this script trains a `MEGNet` model taking `pymatgen` `Structures` as input.

* `explore-structure.py` - this script was intended to plot the bandgap as a function of composition for enumerated disordered variants in the Si-Ge chemical system under the ssumption that we have a solid solution that maintains the ground state Si structure. This is implemented using the`EnumerateStructureTransformation` from `pymatgen`. We produce two plots one where the structures are unscaled and another where we scale based on the estimated volume.

## Scientific Discussion

We cleaned our queried data using the following criteria:
1. We only queried entries for which high quality bandstructure calculations had been carried out.
2. We then discarded all bar the lowest energy polymorph for each composition - this ensures a reasonably defined map for compositon-based models
3. We then discarded any structures that had isolated atoms according to a 4 \AA cutoff criteria.

We performed a stratefied 80:20 train-test split first separating out metals and non-metals and splitting separately before recombining.
We trained two Composition based Random Forest models, one on both metals and non-metals and another on just non-metals.
If we compare how well these two models perform on the non-metals we see that they offer similar performance.


|          | Composition Non-metals |    Composition All    |
|----------|------------------------|-----------------------|
|   MAE    |         0.5026         |        0.5293         |
|   RMSE   |         0.7237         |        0.7422         |
|   R2     |         0.7844         |        0.7733         |

However, these metrics hide the fact that compositon based approaches are ill-suited to the task of bandgap prediction.
This can be seen by looking at the two figures "SiGe-composition-combined.pdf" and "SiGe-composition-non-metals.pdf".
Here we see how ML models are at their heart interpolation devices and the inclusion of the metalic data in the training set leads
to very different performance - note that most of the Si-Ge systems are present in the training set which will have an impact here.

In principle structure-based models should be able to overcome this issue as we have structural degrees of freedom that should
allow the model to accomodate metallic and insulating polymorphs provided they have the sufficiently distinct structures.
Whilst we provided a script capable of training a `MEGNet` model we did not train one in this instance due to available time and resources (Whilst google colab is good for prototyping it is not well suited to long training runs due to the fact that instances are terminate when not interacted with for long periods of time i.e. training a deep model from scratch.)

Due to computational resources I have opted to use the pretrained band gap model provided by `MEGNet` for futher exploration of how structure based ML models perform in the Si-Ge chemical system.
The original authors claim this model has a MAE of 0.33 for their test split of the Materials Project data base.
The Si-Ge chemcial system is as both Si and Ge adopt the same ground state structure and there are no other stable structures in the phase diagram at the DFT level using a
PBE functional.
The experimental band gap (as measured and parameterised in [2]) shows a cross over in behaviour at Si0.15Ge0.85.
We want to test whether the limited data in the Si-Ge system available in the Materials Project is sufficient to recover this beahviour.
The expectation is that the structure-based model will fail to reflect the experimental reality as it is trained on a sparse set of DFT calculated band gaps
that do not reflect experiment themselves.
Most notably the ground state for Ge in the Materials Project is calculated to be metallic.

As structure-based neural networks are extensible and pool over local environments we should in principle be able to estimate a range of possible bandgaps at different compositions by predicting the band gap of different random ordered-disordered superstructures.
I attempted to implement this investigation by using the `EnumerateStructureTransformation` implemented in `pymatgen` however I encountered errors whilst compiling the `enumlib` package used by the transformation class.
I managed to resolve this issue by compiling manually and have written to the authors about updating the conda-forge package to update their recipe. As the fix required manually tweaking several things this is not reproducible but should be once the conda-forge release is updated.

The results from this section are shown in the two figures "SiGe-disordered-structures.pdf" and "SiGe-disordered-structures-scaled.pdf".
The first figure shows the predictions we get without scaling the volume showing a V-shaped plot which badly disagrees with the experimental fit. The second figure shows looks at how scaling the lattice to be approximately the correct volume changes the performance. We see that it infact leads to worse performance.

The conclusion from the rewote exercise is that the data available in the Materials Project is not sufficient to meet the expectation that "the code shall be able to suggest realistic values for slightly modified geometry sets - eg. trained on Si and Ge it should suggest the value of bandgap for Si49Ge51 to be between those of Si and Ge".


## Directions for Extension

1. From a completeness perspective this code could be improved by carrying out hyperparameter searching. However, `RandomForestRegressors` are generally relatively robust to hyper-parameter selection. Consequently, the added value from this in this case is low due to the fact that the composition-based model is ill-suited to the task. Alternatively packages such as `automatminer` which does autoML based on the `tpot` autoML package should be used. These are robust and open sourced
2. Likewise we could in-principle implement hyperparameter searching for the `MEGNet` model but as training deep neural networks is much more expensive and time consuming there was not really a sufficient value proposition for hyperparameter seaching at this time. Accross the literature `MPNN` such as `MEGNet` tends to be relatively robust to selection of hyperparameters - this can be seen from only the slight differences in performance between the 2018 and 2019 pre-trained models that have different hyperparameters.
3. Whilst `MEGNet` satisfies the requirement of being easily adaptable to other materials properties it is not well-suited for the calculation of forces. This is because `MEGNet` uses distances calculated in `pymatgen` outside of the `tensorflow` computation graph. In general if the aim was to build a ML-PES/FF then other models are better suited to the task. Importantly many alternatives already have LAMMPS interfaces making their adoption much easier - Notable examples include Gaussian Approximation Potentials based of the SOAP descriptor, `Flare++` based on using a linearised GP on Atomic Cluster Expansion features, or `NEQUIP` which is an equivariant `MPNN`.
4. `MEGNet` takes a long time for each training epoch (~13 minutes) as the neighbour lists are evaluated through the `GraphConverter` each time an item of data is called. Alternative models opt to pre-computed neighbour lists or cache them to speed up training.
5. Provide useage instructions to allow the work to be reproduced once the update to conda-forge is made.

## Time breakdown

* 2 hours were spent reviewing literature material for Si-Ge and choosing which Machine Learning models to use.
* 2 hours were spent exploring test models and ideas in a google colab notebook before the indermediate meeting.
* 2 hours were spent writing/transfering these ideas into the scripts provided in this repo.
* 1 and a half hours were spent trying to debug an the `enumlib` issue that was blocking the planned deliverables
* 1 and a half hours were spent writing the dicussion and tweaking the scripts.

Total: 9 hrs

## References

[1] https://www.nature.com/articles/npjcompumats201628
[2] https://journals.aps.org/prb/abstract/10.1103/PhysRevB.40.5683
