import math


class NormalizationEncoder(object):
    """Stores the normalization information for a feature.

    The 'factor' data member stores the ratio of the raw variance to the
    normalized variance. The 'encode' method normalizes a value from the
    initial distribution. The 'decode' method de-normalizes a normalized
    value back to the initial distribution.
    """

    def __init__(self, minimum_value: float, maximum_value: float):
        """Create a new NormalizationEncoder instance.

        The first parameter should be the smallest numerical value for a
        given feature in the training data. The second parameter should
        be the largest numerical value for a given feature in the
        training data.
        """
        try:
            minimum_value = float(minimum_value)
        except:
            raise TypeError(
                'Expected type that could be cast to float. Found: ' + str(type(minimum_value)))
        try:
            maximum_value = float(maximum_value)
        except:
            raise TypeError(
                'Expected type that could be case to float. Found: ' + str(type(maximum_value)))
        if not math.isfinite(minimum_value):
            raise ValueError(
                'Expected finite value. Found: ' + str(minimum_value))
        if not math.isfinite(maximum_value):
            raise ValueError(
                'Expected finite value. Found: ' + str(maximum_value))
        if minimum_value == maximum_value:
            self.factor = 1.0
            self.encode = lambda x: 0
            self.decode = lambda x: minimum_value
        elif minimum_value < maximum_value:
            delta = maximum_value - minimum_value
            self.factor = 1.0 / delta
            self.encode = lambda x: (x - minimum_value) / delta
            self.decode = lambda x: x * delta + minimum_value
        else:
            raise ValueError(
                'Expected minimum value to be less than or equal to maximum value. Found: ' + str((minimum_value, maximum_value)))
