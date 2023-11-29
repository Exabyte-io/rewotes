from pymatgen.core import Structure, Composition
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from typing import Union, List
from .featurizer import CompositionFeaturizer, StructureFeaturizer


class Predictor(Pipeline):
    """
    Materials property predictor that predict target property from input structures or stoichiometry.
    """

    def __init__(
            self,
            featurizer: Union[CompositionFeaturizer, StructureFeaturizer] = CompositionFeaturizer(),
            model=LinearRegression(),
            question_type="regression",
    ):
        self.featurizer = featurizer
        self.model = model
        self.question_type = question_type
        steps = [
            (i.__class__.__name__, i)
            for i in [
                self.featurizer,
                self.model,
            ]
            if i
        ]
        super().__init__(steps)