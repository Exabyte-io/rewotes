import hashlib
import numpy
import torch
from ..material import MaterialArchive
from .DataRangeEncoderArray import DataRangeEncoderArray

# int -> bool
# pure function which returns True for 20% of input values
def assign_as_testing(index):
    result = hashlib.md5(str(index).encode("utf-8")).digest()
    return result[0] < 241 and result[1] < 119 and result[2] < 117

def separate_y_x(data):
    y = torch.from_numpy(data[:,0])
    x = torch.from_numpy(data[:,1:])
    return y, x

class Dataset(object):
    def __init__(self, archive):
        if not isinstance(archive, MaterialArchive):
            raise TypeError("Expected instance of MaterialArchive. Found: " + str(type(archive)))
        if len(archive) < 3:
            raise ValueError("Must provide at least three entries for meaningful machine learning. Found: " + str(len(archive)))
        numpy_double_array = archive.to_numpy_double_array()
        self.data_range_encoder_array = DataRangeEncoderArray(numpy_double_array)
        training_data = []
        testing_data = []
        for material_index in range(len(numpy_double_array)):
            material = self.data_range_encoder_array.encode(numpy_double_array[material_index])
            [training_data,testing_data][assign_as_testing(material_index)].append(material)
        self.training_y, self.training_x = separate_y_x(numpy.array(training_data))
        self.testing_y, self.testing_x = separate_y_x(numpy.array(testing_data))
    @staticmethod
    def load_from_archive_file(filename):
        return Dataset(MaterialArchive.load_from_file(filename))
