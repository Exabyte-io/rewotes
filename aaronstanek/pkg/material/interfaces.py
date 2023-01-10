from __future__ import annotations
from abc import ABC
import numpy


class MaterialInterface(ABC):
    """Stores information about a material."""

    def serialize(self) -> bytes:
        """Convert the object to a binary format."""
        raise NotImplementedError()

    def save_to_file(self, filename: str) -> None:
        """Save this object to a file in its binary format."""
        raise NotImplementedError()

    @staticmethod
    def load_from_file(filename: str) -> MaterialInterface:
        """Load this object's binary format from a file, and reconstruct the
        object in memory."""
        raise NotImplementedError()

    def to_numpy_double_array(self) -> numpy.ndarray:
        """Convert this object to 1D numpy array of doubles."""
        raise NotImplementedError()


class MaterialArchiveIteratorInterface(ABC):
    """Iterator type for a given implementation of MaterialArchiveInterface."""

    def __iter__(self) -> MaterialArchiveIteratorInterface:
        raise NotImplementedError()

    def __next__(self) -> MaterialInterface:
        raise NotImplementedError()


class MaterialArchiveInterface(ABC):
    """Linear collection of objects implementing the MaterialInterface."""

    def __len__(self) -> int:
        """Return the number of materials in the collection."""
        raise NotImplementedError()

    def __iter__(self) -> MaterialArchiveIteratorInterface:
        """Return an iterator for this collection."""
        raise NotImplementedError()

    def serialize(self) -> bytes:
        """Convert this collection to a binary format."""
        raise NotImplementedError()

    def __getitem__(self, index: int) -> MaterialInterface:
        """Get the nth material in the collection."""
        raise NotImplementedError()

    def save_to_file(self, filename: str) -> None:
        """Save this collection to a file in its binary format."""
        raise NotImplementedError()

    @staticmethod
    def load_from_file(filename) -> MaterialArchiveInterface:
        """Load this collection's binary format from a file, reconstructing the
        collection in memory."""
        raise NotImplementedError()

    def to_numpy_double_array(self) -> numpy.ndarray:
        """Convert this collection to a 2D numpy array, where each row
        represents the properties of a material in the collection."""
        raise NotImplementedError()

    def save_as_numpy(self, filename: str) -> None:
        """Convert this collection to a 2D numpy array and save the result to a
        file."""
        raise NotImplementedError()
