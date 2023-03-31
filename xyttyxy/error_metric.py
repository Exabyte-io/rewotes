from abc import ABC, abstractmethod


class ErrorMetric(ABC):
    """Abstract class for error metric"""

    def __init__(self, **kwargs):
        if "fractional" in kwargs.keys():
            self.fractional = kwargs["fractional"]
        else:
            self.fractional = False

    @property
    @abstractmethod
    def error(self, calc_0, calc_1):
        raise


class ErrorMetricScalar(ErrorMetric):
    """Error metric using a scalar number, e.g. total energy"""

    def __init__(self, **kwargs):
        ErrorMetric.__init__(self, **kwargs)

    def error(self, calc_0, calc_1):
        v0 = calc_0.raw_value
        v1 = calc_1.raw_value
        v_diff = abs(v0 - v1)

        if self.fractional:
            return v_diff / abs(v0)
        else:
            return v_diff


class ErrorMetricVector(ErrorMetric):
    """Error metric using a vector, e.g. phonon frequencies"""

    def __init__(self, **kwargs):
        ErrorMetric.__init__(self, **kwargs)

    def error(self, calc_0, calc_1):
        raise NotImplementedError


class ErrorMetricMatrix(ErrorMetric):
    """Error metric using a matrix, e.g. x,y,z forces on atoms"""

    def __init__(self, **kwargs):
        ErrorMetric.__init__(self, **kwargs)

    def error(self, calc_0, calc_1):
        raise NotImplementedError
