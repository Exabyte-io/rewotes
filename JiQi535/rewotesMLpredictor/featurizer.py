import numpy as np

from maml.describers._m3gnet import BaseDescriber, M3GNetStructure
from matminer.featurizers.composition import ElementProperty
from pymatgen.core import Structure, Composition
from sklearn.base import BaseEstimator, TransformerMixin
from typing import Union, List

MATMINER_ELEMENT_PRESETS = ["magpie", "deml", "matminer", "matscholar_el", "megnet_el"]


class CompositionFeaturizer(BaseEstimator, TransformerMixin):
    """
    Featurizer encodes composition into features.
    """

    def __init__(self, featurizers: List[ElementProperty] = []):
        """
        Initialize a CompositionFeaturizer.
        :param featurizers: A list of compositional featurizers. Only allows ElementProperty in matminer
        (https://github.com/hackingmaterials/matminer) for the current implementation. Default is Element
        Property with "megnet_el" as preset.
        """
        if featurizers:
            if not (isinstance(featurizers, list) and
                    all(type(x) == ElementProperty for x in featurizers)):
                raise TypeError(
                    f"Featurizers of CompositionFeaturizer must be a list of matminer ElementProperty object."
                    "Check out matminer: https://github.com/hackingmaterials/matminer."
                )
        else:
            featurizers = [ElementProperty.from_preset("megnet_el")]
        self.featurizers = featurizers
        self.feature_labels = [l for f in self.featurizers for l in f.feature_labels()]

    @classmethod
    def from_presets(cls, presets: List[str] = ["megnet_el"]):
        """
        Initialize a CompositionFeaturizer from a list of presets in matminer.
        :param presets: A list as a subset of ["magpie", "deml", "matminer", "matscholar_el", "megnet_el"].
        Default to ["megnet_el"].
        :return: CompositionFeaturizer
        """
        if not presets:
            raise ValueError(
                f"At least one of the presets should be provided."
                f"Options: {MATMINER_ELEMENT_PRESETS}."
            )
        for preset in presets:
            if preset not in MATMINER_ELEMENT_PRESETS:
                raise ValueError(f"{preset} is not an allowed preset in matminer.")
        featurizers = [ElementProperty.from_preset(p) for p in presets]
        return CompositionFeaturizer(featurizers=featurizers)

    def fit(self, X = None, y = None):
        """
        Place holder for fit method.
        :param X: Any input.
        :param y: Any input.
        :return: self
        """
        return self

    def transform(self, compositions: Union[List[str], List[Composition]] = None):
        """
        Transform compositions into features.
        :param compositions: A list of str or a list of pymatgen Composition.
        :return: A numpy array of compositional features.
        """
        compositions = self._check_compositions(compositions)
        composition_features = np.concatenate([f.transform(compositions) for f in self.featurizers], axis=1)
        return composition_features

    def predict(self, compositions: Union[List[str], List[Composition]] = None):
        """
        Transform compositions into features. Same as transform method.
        :param compositions: A list of str or a list of pymatgen Composition.
        :return: A numpy array of compositional features.
        """
        return self.transform(compositions)

    @staticmethod
    def _check_compositions(compositions: Union[List[str], List[Composition]] = None):
        """
        Check the compositions are in correct format.
        :param compositions: A list of str or a list of pymatgen Composition.
        :return: A list of pymatgen Composition.
        """
        if not (isinstance(compositions, list) and
                all(isinstance(x, (str, Composition)) for x in compositions)):
            raise TypeError("Compositions must be provided as a list of str or a list of pymatgen Composition.")
        if any(type(x) == str for x in compositions):
            compositions = [Composition(str(x)) for x in compositions]
        return compositions


class StructureFeaturizer(BaseEstimator, TransformerMixin):
    """
    Featurizer encodes structure into features.
    """

    def __init__(self, featurizers: List[BaseDescriber] = []):
        """
        Initialize a StructureFeaturizer.
        :param featurizers: A list of structural featurizers. Only allows BaseDescriber in maml
        (https://github.com/materialsvirtuallab/maml) for the current implementation. Default is
        M3GNet structure featurizer.
        """
        if featurizers:
            if not (isinstance(featurizers, list) and
                    all(type(x) == BaseDescriber for x in featurizers)):
                raise TypeError(
                    f"Featurizers of StructureFeaturizer must be a list of maml BaseDescriber object."
                    "Check out maml: https://github.com/materialsvirtuallab/maml."
                )
        else:
            featurizers = [M3GNetStructure()]
        self.featurizers = featurizers

    def fit(self, X = None, y = None):
        """
        Place holder for fit method.
        :param X: Any input.
        :param y: Any input.
        :return: self
        """
        return self

    def transform(self, structures: List[Structure] = None):
        """
        Transform structures into features.
        :param structures: A list of pymatgen structures.
        :return: A numpy array of structural features.
        """
        structures = self._check_structures(structures)
        features = np.concatenate([f.transform(structures) for f in self.featurizers], axis=1)
        return features

    def predict(self, structures: List[Structure] = None):
        """
        Transform structures into features. Same as transform method.
        :param structures: A list of pymatgen structures.
        :return:  A numpy array of structural features.
        """
        return self.transform(structures)

    @staticmethod
    def _check_structures(structures: List[Structure] = None):
        """
        Check the structures are in correct format.
        :param structures: A list of pymatgen Composition.
        :return: A list of pymatgen structures.
        """
        if not (isinstance(structures, list) and
                all(type(x) == Structure for x in structures)):
            raise TypeError("Structures must be provided as a list of pymatgen Structures.")
        return structures
