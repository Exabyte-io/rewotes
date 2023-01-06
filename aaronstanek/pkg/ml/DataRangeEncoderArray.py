import numpy
from .DataRangeEncoder import DataRangeEncoder

class DataRangeEncoderArray(object):
    def __init__(self, numpy_double_array):
        if type(numpy_double_array) != numpy.ndarray:
            raise TypeError("Expected numpy array. Found: " + str(type(numpy_double_array)))
        if len(numpy_double_array.shape) != 2:
            raise ValueError("Expected numpy array to be 2D. Found shape: " + str(numpy_double_array.shape))
        minimum_values = numpy.apply_along_axis(min, 1, numpy_double_array)
        maximum_values = numpy.apply_along_axis(max, 1, numpy_double_array)
        self.encoders = list(map(
            lambda property_index: DataRangeEncoder(minimum_values[property_index], maximum_values[property_index]),
            range(numpy_double_array.shape[1])
        ))
    def __len__(self):
        return len(self.encoders)
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
