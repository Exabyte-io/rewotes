class metrix():
    """Measures error of convergence testing"""

    def __init__(self):
        pass

    def __call__(self, value_1: float, value_2: float) -> float:
        """Measures error between two points"""
        pass


class mae(metrix):
    """Mean absolute error for convergence testing"""

    def __call__(self, value_1: float, value_2: float) -> float:
        """Measures error between two points as mean absolute error"""
        return abs(value_1 - value_2)