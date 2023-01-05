import numpy
from .material_pb2 import MaterialArchive as protobuf_material_archive
from .Material import Material

class MaterialArchiveIterator(object):
    def __init__(self, archive):
        self.archive = archive
        self.next_index = 0
    def __iter__(self):
        return self
    def __next__(self):
        if self.next_index >= len(self.archive):
            raise StopIteration
        else:
            return_value = self.archive[self.next_index]
            self.next_index += 1
            return return_value

class MaterialArchive(object):
    def __init__(self, initializer = None):
        self._material_archive = protobuf_material_archive()
        if initializer is None:
            pass
        elif type(initializer) == bytes:
            self._material_archive.ParseFromString(initializer)
        else:
            raise TypeError("Expected bytes. Found: " + str(type(initializer)))
    def __len__(self):
        return len(self._material_archive.serialized_materials)
    def __iter__(self):
        return MaterialArchiveIterator(self)
    def __getitem__(self, index):
        try:
            index = int(index)
        except:
            raise TypeError("Expected type that could be cast to int. Found: " + str(type(index)))
        return Material(self._material_archive.serialized_materials[index])
    def append(self, material):
        if not isinstance(material, Material):
            raise TypeError("Expected instance of Material. Found: " + str(type(material)))
        self._material_archive.serialized_materials.append(material.serialize())
    def serialize(self):
        return self._material_archive.SerializeToString()
    def save_to_file(self, filename):
        if type(filename) != str:
            raise TypeError("Expected str. Found: " + str(type(filename)))
        else:
            with open(filename, "wb") as file:
                file.write(self.serialize())
    @staticmethod
    def load_from_file(filename):
        if type(filename) != str:
            raise TypeError("Expected str. Found: " + str(type(filename)))
        else:
            with open(filename, "rb") as file:
                return MaterialArchive(file.read())
    def to_numpy_double_array(self):
        rows = []
        for material in self:
            rows.append(material.to_numpy_double_array())
        return numpy.array(rows)
