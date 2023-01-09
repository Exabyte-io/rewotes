from __future__ import annotations
import numpy
from .material_pb2 import Material as protobuf_material
from typing import Any, Optional
from .interfaces import MaterialInterface


class Material(MaterialInterface):
    """Thin wrapper around Material protocol-buffer."""

    def __init__(self, initializer: Optional[bytes] = None):
        """Create a new Material instance.

        Optionally pass a serialized Material protocol-buffer as
        initializer parameter.
        """
        self.__dict__['_values'] = protobuf_material()
        if initializer is None:
            pass
        elif type(initializer) == bytes:
            self.__dict__['_values'].ParseFromString(initializer)
        else:
            raise TypeError('Expected bytes. Found: ' + str(type(initializer)))

    def __getattr__(self, property_name: str) -> Any:
        """Get the value of a property stored in the wrapped protocol-
        buffer."""
        return getattr(self.__dict__['_values'], property_name)

    def __setattr__(self, property_name: str, property_value: Any) -> None:
        """Set the value of a property in the wrapped protocol-buffer."""
        setattr(self.__dict__['_values'], property_name, property_value)

    def __delattr__(self, property_name: str) -> None:
        """Attempt to delete a property in the wrapped protocol-buffer."""
        delattr(self.__dict__['_values'], property_name)

    def serialize(self) -> bytes:
        """Return a serialized version the wrapped protocol-buffer."""
        return self.__dict__['_values'] .SerializeToString()

    def save_to_file(self, filename: str) -> None:
        """Serialize the wrapped protocol-buffer and save the result to a
        file."""
        if type(filename) != str:
            raise TypeError('Expected str. Found: ' + str(type(filename)))
        else:
            with open(filename, 'wb') as file:
                file.write(self.serialize())

    @staticmethod
    def load_from_file(filename: str) -> Material:
        """Read a file, interpret the contents as a serialized Material
        protocol-buffer, and return an instance of Material wrapping the de-
        serialized protocol-buffer."""
        if type(filename) != str:
            raise TypeError('Expected str. Found: ' + str(type(filename)))
        else:
            with open(filename, 'rb') as file:
                return Material(file.read())

    def to_numpy_double_array(self) -> numpy.ndarray:
        """Create a 1D numpy array with dtype=double from the entire contents
        of the wrapped protocol-buffer.

        THE ADDITION OR REMOVAL OF PROPERTIES IN THE MATERIAL PROTOCOL-BUFFER WILL NECESSITATE BREAKING CHANGES IN THIS METHOD!
        Within any given version of this module, the arrays returned from this method will have consistent length.
        Within any given version of this module, the interpretation of the value at each index will be consistent.
        The value at index = 0 will always represent the bang gap of the material.
        """
        values = []
        values.append(self.band_gap)
        for atomic_number in range(1, 119):
            values.append(self.composition[atomic_number])
        for atomic_number in range(1, 119):
            values.append(self.composition_reduced[atomic_number])
        return numpy.array(values, dtype=numpy.double)
