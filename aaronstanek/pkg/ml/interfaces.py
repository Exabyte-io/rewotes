from __future__ import annotations
from abc import ABC
from typing import Any, Tuple
import numpy


class DatasetInterface(ABC):
    '''Collection of feature rows with their corresponding target values.'''
    def __len__(self) -> int:
        '''Return the number of rows in the collection.'''
        raise NotImplementedError()

    def __getitem__(self, index: int) -> Tuple[Any, Any]:
        '''Return the nth feature row and the nth target value as a tuple.'''
        raise NotImplementedError()


class NormalizationEncoderInterface(ABC):
    '''A type to store the conversion functions between the normalized and unnormalized form of a feature.'''
    def get_factor(self) -> float:
        '''Return the result of the raw variance divided by the normalized variance.'''
        raise NotImplementedError()

    def encode(self, x: float) -> float:
        '''Normalize a value.'''
        raise NotImplementedError()

    def decode(self, x: float) -> float:
        '''De-normalize a value.'''
        raise NotImplementedError()


class NormalizationEncoderArrayInterface(ABC):
    '''Linear collection of objects implementing the NormalizationEncoderInterface.'''
    def __len__(self) -> int:
        '''Return the number of encoders in the collection.'''
        raise NotImplementedError()

    def __getitem__(self, index: int) -> NormalizationEncoderInterface:
        '''Return the nth encoder.'''
        raise NotImplementedError()

    def encode(self, numpy_double_array: numpy.ndarray) -> numpy.ndarray:
        '''Normalize a feature row, where the nth value is normalized by the nth encoder.'''
        raise NotImplementedError()

    def decode(self, numpy_double_array: numpy.ndarray) -> numpy.ndarray:
        '''De-normalize a feature row, where the nth value is de-normalized by the nth encoder.'''
        raise NotImplementedError()


class TrainingDataManagerInterface(ABC):
    """Builds and holds objects that implement the Dataset Interface objects."""
    @staticmethod
    def load_from_archive_file(filename: str) -> TrainingDataManagerInterface:
        '''Load a binary file created by a type implementing MaterialArchiveInterface.'''
        raise NotImplementedError()

    @staticmethod
    def load_from_numpy_file(filename: str) -> TrainingDataManagerInterface:
        '''Load a numpy file created by a type implementing MaterialArchiveInterface.'''
        raise NotImplementedError()
