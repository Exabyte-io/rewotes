from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from typing import Union
from .featurizer import CompositionFeaturizer, StructureFeaturizer


class Predictor(Pipeline):
    """
    Materials property predictor that predict target property from input structures or stoichiometry.
    """

    def __init__(
            self,
            featurizer: Union[CompositionFeaturizer, StructureFeaturizer] = CompositionFeaturizer(),
            model=LinearRegression(),
    ):
        """
        Initialize the ML property predictors.
        :param featurizer: CompositionFeaturizer or StructureFeaturizer. Default is CompositionFeaturizer.
        :param model: The regression model in sklearn to use. Default is LinearRegression model.
        """
        self.featurizer = featurizer
        self.model = model
        steps = [
            (i.__class__.__name__, i)
            for i in [
                self.featurizer,
                self.model,
            ]
            if i
        ]
        super().__init__(steps)