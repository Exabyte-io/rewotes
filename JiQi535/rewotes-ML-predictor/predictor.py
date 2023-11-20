from sklearn.pipeline import Pipeline

class Predictor(Pipeline):
    """
    Materials property predictor that predict target property from input structures or stoichiometry.
    """
    def __init__(self, model=None):
        self.model=model

    def fit(self, X=None, y=None):
        return False

    def transform(self, X):
        return X