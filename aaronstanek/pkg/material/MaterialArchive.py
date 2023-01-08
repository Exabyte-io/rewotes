from __future__ import annotations
import numpy
from .material_pb2 import MaterialArchive as protobuf_material_archive
from .Material import Material
from typing import Optional

class MaterialArchiveIterator(object):
    def __init__(self, archive):
        self.archive = archive
        self.next_index = 0
    def __iter__(self) -> MaterialArchiveIterator:
        return self
    def __next__(self) -> Material:
        if self.next_index >= len(self.archive):
            raise StopIteration
        else:
            return_value = self.archive[self.next_index]
            self.next_index += 1
            return return_value

class MaterialArchive(object):
    def __init__(self, initializer: Optional[bytes] = None):
        self._material_archive = protobuf_material_archive()
        if initializer is None:
            pass
        elif type(initializer) == bytes:
            self._material_archive.ParseFromString(initializer)
        else:
            raise TypeError("Expected bytes. Found: " + str(type(initializer)))
    def __len__(self) -> int:
        return len(self._material_archive.serialized_materials)
    def __iter__(self) -> MaterialArchiveIterator:
        return MaterialArchiveIterator(self)
    def __getitem__(self, index: int) -> Material:
        try:
            index = int(index)
        except:
            raise TypeError("Expected type that could be cast to int. Found: " + str(type(index)))
        return Material(self._material_archive.serialized_materials[index])
    def append(self, material: Material) -> None:
        if not isinstance(material, Material):
            raise TypeError("Expected instance of Material. Found: " + str(type(material)))
        self._material_archive.serialized_materials.append(material.serialize())
    def serialize(self) -> bytes:
        return self._material_archive.SerializeToString()
    def save_to_file(self, filename: str) -> None:
        if type(filename) != str:
            raise TypeError("Expected str. Found: " + str(type(filename)))
        else:
            with open(filename, "wb") as file:
                file.write(self.serialize())
    @staticmethod
    def load_from_file(filename) -> MaterialArchive:
        if type(filename) != str:
            raise TypeError("Expected str. Found: " + str(type(filename)))
        else:
            with open(filename, "rb") as file:
                return MaterialArchive(file.read())
    def to_numpy_double_array(self) -> numpy.ndarray:
        rows = []
        for material in self:
            rows.append(material.to_numpy_double_array())
        return numpy.array(rows)
    def save_as_numpy(self, filename: str) -> None:
        if type(filename) != str:
            raise TypeError("Expected str. Found: " + str(type(filename)))
        else:
            numpy_double_array = self.to_numpy_double_array()
            numpy.save(filename, numpy_double_array, allow_pickle=False, fix_imports=False)