import math

class DataRangeEncoder(object):
    def __init__(self, minimum_value, maximum_value):
        try:
            minimum_value = float(minimum_value)
        except:
            raise TypeError("Expected type that could be cast to float. Found: " + str(type(minimum_value)))
        try:
            maximum_value = float(maximum_value)
        except:
            raise TypeError("Expected type that could be case to float. Found: " + str(type(maximum_value)))
        if not math.isfinite(minimum_value):
            raise ValueError("Expected finite value. Found: " + str(minimum_value))
        if not math.isfinite(maximum_value):
            raise ValueError("Expected finite value. Found: " + str(maximum_value))
        if minimum_value == maximum_value:
            self.encode = lambda x: 0
            self.decode = lambda x: minimum_value
        elif minimum_value < maximum_value:
            delta = maximum_value - minimum_value
            self.encode = lambda x: (x - minimum_value) / delta
            self.decode = lambda x: x * delta + minimum_value
        else:
            raise ValueError("Expected minimum value to be less than or equal to maximum value. Found: " + str((minimum_value, maximum_value)))
