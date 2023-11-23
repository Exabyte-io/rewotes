import numpy as np

from maml.base import BaseDescriber
from matminer.featurizers.composition import ElementProperty
from pymatgen.core import Structure, Composition
from sklearn.pipeline import Pipeline
from typing import Union, List

MATMINER_ELEMENT_PRESETS = ["magpie", "deml", "matminer", "matscholar_el", "megnet_el"]


class Featurizer(Pipeline):
    """
    Featurizer encodes structure or stoichiometry into features.
    """

    def __init__(self, composition_featurizers: List[ElementProperty] = [], structure_featurizers: list = []):
        if composition_featurizers:
            if not (isinstance(composition_featurizers, list) and
                    all(type(x) == ElementProperty for x in composition_featurizers)):
                raise TypeError(
                    f"composition_featurizers must be a list of matminer ElementProperty object."
                    "Check out matminer: https://github.com/hackingmaterials/matminer."
                )
        self.composition_featurizers = composition_featurizers
        self.feature_labels = [l for f in self.composition_featurizers for l in f.feature_labels()]

        if structure_featurizers:
            if not (isinstance(structure_featurizers, list) and
                    all(type(x) == BaseDescriber for x in structure_featurizers)):
                raise TypeError(
                    f"structure_featurizers must be a list of maml BaseDescriber object."
                    "Check out maml: https://github.com/materialsvirtuallab/maml."
                )
        self.structure_featurizers = structure_featurizers
        # feature labels to be defined.
        # self.feature_labels =

    def fit(self, compositions: Union[List[str], List[Composition]] = None, structures: List[Structure] = None):
        return self.transform(compositions, structures)

    def transform(self, compositions: Union[List[str], List[Composition]] = None, structures: List[Structure] = None):
        structures, compositions = self._check_structure_and_compositions(compositions, structures)
        if self.composition_featurizers:
            composition_features = np.concatenate([f.transform(compositions) for f in self.composition_featurizers],
                                                  axis=1)
            return composition_features

        # structure feature transform to be defined

    def _check_structure_and_compositions(self, compositions: Union[List[str], List[Composition]] = None,
                                          structures: List[Structure] = None):
        if structures:
            if not (isinstance(structures, list) and
                    all(type(x) == Structure for x in structures)):
                raise TypeError("Structures must be provided as a list of pymatgen Structures.")
        elif self.structure_featurizers:
            raise ValueError("To process structure_features, structures must be provided.")

        if compositions:
            if not (isinstance(compositions, list) and
                    all(isinstance(x, (str, Composition)) for x in compositions)):
                raise TypeError("Compositions must be provided as a list of str or a list of pymatgen Composition.")
            if any(type(x) == str for x in compositions):
                compositions = [Composition(str(x)) for x in compositions]
        elif self.composition_featurizers:
            compositions = [s.composition for s in structures]

        return structures, compositions

class CompositionFeaturizer(Featurizer):
    def __init__(self, presets: List[str] = ['megnet_el']):
        for preset in presets:
            if preset not in MATMINER_ELEMENT_PRESETS:
                raise ValueError(f"{preset} is not an allowed preset in matminer.")
        composition_featurizers = [ElementProperty.from_preset(p) for p in presets]
        super().__init__(composition_featurizers)
        self.presets = presets
