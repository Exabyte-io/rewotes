import numpy
from .DataRangeEncoder import DataRangeEncoder

class DataRangeEncoderArray(object):
    def __init__(self, numpy_double_array):
        if type(numpy_double_array) != numpy.ndarray:
            raise TypeError("Expected numpy array. Found: " + str(type(numpy_double_array)))
        if len(numpy_double_array.shape) != 2:
            raise ValueError("Expected numpy array to be 2D. Found shape: " + str(numpy_double_array.shape))
        minimum_values = []
        maximum_values = []
        for material in numpy_double_array:
            for property_index in range(0, len(material)):
                if len(minimum_values) > property_index:
                    minimum_values[property_index] = min(minimum_values[property_index], material[property_index])
                    maximum_values[property_index] = max(maximum_values[property_index], material[property_index])
                else:
                    minimum_values.append(material[property_index])
                    maximum_values.append(material[property_index])
        self.encoders = []
        for property_index in range(len(minimum_values)):
            self.encoders.append(DataRangeEncoder(minimum_values[property_index], maximum_values[property_index]))
    def encode(self, numpy_double_array):
        if type(numpy_double_array) != numpy.ndarray:
            raise TypeError("Expected numpy array. Found: " + str(type(numpy_double_array)))
        if len(numpy_double_array.shape) != 1:
            raise ValueError("Expected numpy array to be 1D. Found shape: " + str(numpy_double_array.shape))
        if len(numpy_double_array) != len(self.encoders):
            raise ValueError("Expected numpy array to have the same length as the sample. Found: " + str(len(numpy_double_array)) + " != " + str(len(self.encoders)))
        return numpy.array(list(map(
            lambda property_index: self.encoders[property_index].encode(numpy_double_array[property_index]),
            range(len(self.encoders)))),
            dtype=numpy.double)
    def decode(self, numpy_double_array):
        if type(numpy_double_array) != numpy.ndarray:
            raise TypeError("Expected numpy array. Found: " + str(type(numpy_double_array)))
        if len(numpy_double_array.shape) != 1:
            raise ValueError("Expected numpy array to be 1D. Found shape: " + str(numpy_double_array.shape))
        if len(numpy_double_array) != len(self.encoders):
            raise ValueError("Expected numpy array to have the same length as the sample. Found: " + str(len(numpy_double_array)) + " != " + str(len(self.encoders)))
        return numpy.array(list(map(
            lambda property_index: self.encoders[property_index].decode(numpy_double_array[property_index]),
            range(len(self.encoders)))),
            dtype=numpy.double)
