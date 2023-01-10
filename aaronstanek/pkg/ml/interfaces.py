from __future__ import annotations
from abc import ABC
from typing import Any
import numpy


class DatasetInterface(ABC):
    def __len__(self) -> int:
        raise NotImplementedError()

    def __getitem__(self, index: int) -> Any:
        raise NotImplementedError()


class NormalizationEncoderInterface(ABC):
    def get_factor(self) -> float:
        raise NotImplementedError()

    def encode(self, x: float) -> float:
        raise NotImplementedError()

    def decode(self, x: float) -> float:
        raise NotImplementedError()


class NormalizationEncoderArrayInterface(ABC):
    def __len__(self) -> int:
        raise NotImplementedError()

    def __getitem__(self, index: int) -> NormalizationEncoderInterface:
        raise NotImplementedError()

    def encode(self, numpy_double_array: numpy.ndarray) -> numpy.ndarray:
        raise NotImplementedError()

    def decode(self, numpy_double_array: numpy.ndarray) -> numpy.ndarray:
        raise NotImplementedError()


class TrainingDataManagerInterface(ABC):
    @staticmethod
    def load_from_archive_file(filename: str) -> TrainingDataManagerInterface:
        raise NotImplementedError()

    @staticmethod
    def load_from_numpy_file(filename: str) -> TrainingDataManagerInterface:
        raise NotImplementedError()
