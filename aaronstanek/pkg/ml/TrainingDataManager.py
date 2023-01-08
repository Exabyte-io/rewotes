from __future__ import annotations
import hashlib
import numpy
import torch
from ..material import MaterialArchive
from .DataRangeEncoderArray import DataRangeEncoderArray
from .Dataset import Dataset
from typing import Union

# int -> bool
# pure function which returns True for 20% of input values
def assign_as_testing(index: int) -> bool:
    result = hashlib.md5(str(index).encode("utf-8")).digest()
    return result[0] < 241 and result[1] < 119 and result[2] < 117

class TrainingDataManager(object):
    def __init__(self, initializer: Union[MaterialArchive, numpy.ndarray], batch_size: int = 64):
        try:
            batch_size = int(batch_size)
        except:
            raise TypeError("Expected type that could be cast to int. Found: " + str(type(batch_size)))
        if isinstance(initializer, MaterialArchive):
            if len(initializer) < 3:
                raise ValueError("Must provide at least three entries for meaningful machine learning. Found: " + str(len(initializer)))
            numpy_double_array = initializer.to_numpy_double_array()
        elif isinstance(initializer, numpy.ndarray):
            if initializer.dtype != numpy.double:
                raise TypeError("Expected numpy array to have dtype=double. Found: " + str(initializer.dtype))
            if initializer.shape[0] < 3:
                raise ValueError("Must provide at least three entries for meaningful machine learning. Found: " + str(initializer.shape[0]))
            if initializer.shape[1] < 2:
                raise ValueError("Must provide at least one feature on which to train.")
            numpy_double_array = initializer
        self.data_range_encoder_array = DataRangeEncoderArray(numpy_double_array)
        self.total_feature_width = len(self.data_range_encoder_array)
        training_data = []
        testing_data = []
        for material_index in range(len(numpy_double_array)):
            material = self.data_range_encoder_array.encode(numpy_double_array[material_index])
            [training_data,testing_data][assign_as_testing(material_index)].append(material)
        self.training = torch.utils.data.DataLoader(Dataset(training_data), batch_size=batch_size, shuffle=True)
        self.testing = torch.utils.data.DataLoader(Dataset(testing_data), batch_size=batch_size, shuffle=True)
    @staticmethod
    def load_from_archive_file(filename: str) -> TrainingDataManager:
        return TrainingDataManager(MaterialArchive.load_from_file(filename))
    @staticmethod
    def load_from_numpy_file(filename: str) -> TrainingDataManager:
        return TrainingDataManager(numpy.load(filename))
