from __future__ import annotations
import numpy
import torch
from ..material import MaterialArchive
from .NormalizationEncoderArray import NormalizationEncoderArray
from .Dataset import Dataset
from typing import Union


class TrainingDataManager(object):
    """Builds and holds training and testing Dataset objects."""

    def __init__(self, initializer: Union[MaterialArchive, numpy.ndarray], batch_size: int = 64):
        """Create a TrainingDataManager from either a MaterialArchive or a
        numpy array."""
        try:
            batch_size = int(batch_size)
        except:
            raise TypeError(
                'Expected type that could be cast to int. Found: ' + str(type(batch_size)))
        if isinstance(initializer, MaterialArchive):
            if len(initializer) < 3:
                raise ValueError(
                    'Must provide at least three entries for meaningful machine learning. Found: ' + str(len(initializer)))
            numpy_double_array = initializer.to_numpy_double_array()
        elif isinstance(initializer, numpy.ndarray):
            if initializer.dtype != numpy.double:
                raise TypeError(
                    'Expected numpy array to have dtype=double. Found: ' + str(initializer.dtype))
            if initializer.shape[0] < 3:
                raise ValueError(
                    'Must provide at least three entries for meaningful machine learning. Found: ' + str(initializer.shape[0]))
            if initializer.shape[1] < 2:
                raise ValueError(
                    'Must provide at least one feature on which to train.')
            numpy_double_array = initializer
        self.data_range_encoder_array = NormalizationEncoderArray(
            numpy_double_array)
        self.total_feature_width = len(self.data_range_encoder_array)
        training_data = []
        testing_data = []
        rng = numpy.random.Generator(numpy.random.MT19937(6))
        for material_index in range(len(numpy_double_array)):
            material = self.data_range_encoder_array.encode(
                numpy_double_array[material_index])
            [training_data, testing_data][int(
                rng.integers(0, 5) == 0)].append(material)
        self.training = torch.utils.data.DataLoader(
            Dataset(training_data), batch_size=batch_size, shuffle=True)
        self.testing = torch.utils.data.DataLoader(
            Dataset(testing_data), batch_size=batch_size, shuffle=True)

    @staticmethod
    def load_from_archive_file(filename: str) -> TrainingDataManager:
        """Create a TrainingDataManager using a file generated from
        MaterialArchive.save_to_file."""
        return TrainingDataManager(MaterialArchive.load_from_file(filename))

    @staticmethod
    def load_from_numpy_file(filename: str) -> TrainingDataManager:
        """Create a TrainingDataManager using a file generated from
        MaterialArchive.save_as_numpy."""
        return TrainingDataManager(numpy.load(filename))
