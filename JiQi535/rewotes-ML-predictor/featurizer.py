from sklearn.pipeline import Pipeline

class Featurizer(Pipeline):
    """
    Featurizer encodes structure or stoichiometry into features.
    """
    def __init__(self, composition_model=True):
        self.composition_model = composition_model

    def fit(self, X=None, y=None):
        return X

    def transform(self, X=None):
        return X