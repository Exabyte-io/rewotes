from __future__ import annotations
from abc import ABC
import numpy


class MaterialInterface(ABC):
    def serialize(self) -> bytes:
        raise NotImplementedError()

    def save_to_file(self, filename: str) -> None:
        raise NotImplementedError()

    @staticmethod
    def load_from_file(filename: str) -> MaterialInterface:
        raise NotImplementedError()

    def to_numpy_double_array(self) -> numpy.ndarray:
        raise NotImplementedError()


class MaterialArchiveIteratorInterface(ABC):
    def __iter__(self) -> MaterialArchiveIteratorInterface:
        raise NotImplementedError()

    def __next__(self) -> MaterialInterface:
        raise NotImplementedError()


class MaterialArchiveInterface(ABC):
    def __len__(self) -> int:
        raise NotImplementedError()

    def __iter__(self) -> MaterialArchiveIteratorInterface:
        raise NotImplementedError()

    def serialize(self) -> bytes:
        raise NotImplementedError()

    def __getitem__(self, index: int) -> MaterialInterface:
        raise NotImplementedError()

    def save_to_file(self, filename: str) -> None:
        raise NotImplementedError()

    @staticmethod
    def load_from_file(filename) -> MaterialArchiveInterface:
        raise NotImplementedError()

    def to_numpy_double_array(self) -> numpy.ndarray:
        raise NotImplementedError()

    def save_as_numpy(self, filename: str) -> None:
        raise NotImplementedError()
