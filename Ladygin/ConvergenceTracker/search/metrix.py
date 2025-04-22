from abc import ABC, abstractmethod


class metrix(ABC):
    """Measures error of convergence testing"""
    
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, value_1, value_2) -> float:
        """Measures error between two points"""
        raise NotImplementedError


class mae(metrix):
    """Mean absolute error for convergence testing"""
    def __init__(self):
        super().__init__()

    
    def __call__(self, value_1, value_2) -> float:
        """Measures error between two points as mean absolute error"""
        try:
            return sum(map(lambda x1, x2: abs(x1 - x2), value_1, value_2))/len(value_1)
        except TypeError:
            return abs(value_1 - value_2)
    