from __future__ import annotations
import numpy
from .material_pb2 import Material as protobuf_material
from typing import Any, Optional


class Material(object):
    def __init__(self, initializer: Optional[bytes] = None):
        self.__dict__['_values'] = protobuf_material()
        if initializer is None:
            pass
        elif type(initializer) == bytes:
            self.__dict__['_values'].ParseFromString(initializer)
        else:
            raise TypeError('Expected bytes. Found: ' + str(type(initializer)))

    def __getattr__(self, property_name: str) -> Any:
        return getattr(self.__dict__['_values'], property_name)

    def __setattr__(self, property_name: str, property_value: Any) -> None:
        setattr(self.__dict__['_values'], property_name, property_value)

    def __delattr__(self, property_name: str) -> None:
        delattr(self.__dict__['_values'], property_name)

    def serialize(self) -> bytes:
        return self.__dict__['_values'] .SerializeToString()

    def save_to_file(self, filename: str) -> None:
        if type(filename) != str:
            raise TypeError('Expected str. Found: ' + str(type(filename)))
        else:
            with open(filename, 'wb') as file:
                file.write(self.serialize())

    @staticmethod
    def load_from_file(filename: str) -> Material:
        if type(filename) != str:
            raise TypeError('Expected str. Found: ' + str(type(filename)))
        else:
            with open(filename, 'rb') as file:
                return Material(file.read())

    def to_numpy_double_array(self) -> numpy.ndarray:
        values = []
        values.append(self.band_gap)
        for atomic_number in range(1, 119):
            values.append(self.composition[atomic_number])
        for atomic_number in range(1, 119):
            values.append(self.composition_reduced[atomic_number])
        return numpy.array(values, dtype=numpy.double)
