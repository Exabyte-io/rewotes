from .interfaces import NormalizationEncoderArrayInterface
import numpy
from .NormalizationEncoder import NormalizationEncoder


class NormalizationEncoderArray(NormalizationEncoderArrayInterface):
    """Collection of NormalizationEncoder objects.

    Provides normalization and de-normalization for all features in a
    dataset.
    """

    def __init__(self, numpy_double_array: numpy.ndarray):
        """Create a new NormalizationEncoderArray object.

        Internally generate and store a NormalizationEncoder for each
        feature in the input numpy array.
        """
        if type(numpy_double_array) != numpy.ndarray:
            raise TypeError('Expected numpy array. Found: ' +
                            str(type(numpy_double_array)))
        if len(numpy_double_array.shape) != 2:
            raise ValueError(
                'Expected numpy array to be 2D. Found shape: ' + str(numpy_double_array.shape))
        minimum_values = numpy.apply_along_axis(min, 0, numpy_double_array)
        maximum_values = numpy.apply_along_axis(max, 0, numpy_double_array)
        self.encoders = list(map(
            lambda property_index: NormalizationEncoder(
                minimum_values[property_index], maximum_values[property_index]),
            range(len(minimum_values))
        ))

    def __len__(self) -> int:
        """Return the number of NormalizationEncoder objects stored in the
        collection.

        This is equal to the sum of the number of input and target
        features.
        """
        return len(self.encoders)

    def __getitem__(self, index: int) -> NormalizationEncoder:
        """Return the encoder for the nth feature."""
        try:
            index = int(index)
        except:
            raise TypeError(
                'Expected type that could be cast to int. Found: ' + str(type(index)))
        return self.encoders[index]

    def encode(self, numpy_double_array: numpy.ndarray) -> numpy.ndarray:
        """Return a normalized copy of an input row."""
        if type(numpy_double_array) != numpy.ndarray:
            raise TypeError('Expected numpy array. Found: ' +
                            str(type(numpy_double_array)))
        if len(numpy_double_array.shape) != 1:
            raise ValueError(
                'Expected numpy array to be 1D. Found shape: ' + str(numpy_double_array.shape))
        if len(numpy_double_array) != len(self.encoders):
            raise ValueError('Expected numpy array to have the same length as the sample. Found: ' +
                             str(len(numpy_double_array)) + ' != ' + str(len(self.encoders)))
        return numpy.array(list(map(
            lambda property_index: self.encoders[property_index].encode(
                numpy_double_array[property_index]),
            range(len(self.encoders)))),
            dtype=numpy.double)

    def decode(self, numpy_double_array: numpy.ndarray) -> numpy.ndarray:
        """Return a de-normalized copy of a normalized input row."""
        if type(numpy_double_array) != numpy.ndarray:
            raise TypeError('Expected numpy array. Found: ' +
                            str(type(numpy_double_array)))
        if len(numpy_double_array.shape) != 1:
            raise ValueError(
                'Expected numpy array to be 1D. Found shape: ' + str(numpy_double_array.shape))
        if len(numpy_double_array) != len(self.encoders):
            raise ValueError('Expected numpy array to have the same length as the sample. Found: ' +
                             str(len(numpy_double_array)) + ' != ' + str(len(self.encoders)))
        return numpy.array(list(map(
            lambda property_index: self.encoders[property_index].decode(
                numpy_double_array[property_index]),
            range(len(self.encoders)))),
            dtype=numpy.double)
