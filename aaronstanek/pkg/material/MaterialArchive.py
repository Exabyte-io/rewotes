from __future__ import annotations
import numpy
from .material_pb2 import MaterialArchive as protobuf_material_archive
from .Material import Material
from typing import Optional
from .interfaces import MaterialArchiveIteratorInterface, MaterialArchiveInterface


class MaterialArchiveIterator(MaterialArchiveIteratorInterface):
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


class MaterialArchive(MaterialArchiveInterface):
    """Memory-efficient collection of Material objects.

    Reading from and mutating this object are slow operations. A Python
    list is faster, but consumes substantially more memory. Internally
    uses a MaterialArchive protocol-buffer.
    """

    def __init__(self, initializer: Optional[bytes] = None):
        """Create a new MaterialArchive instance.

        Optionally pass a serialized MaterialArchive protocol-buffer as
        initializer parameter.
        """
        self._material_archive = protobuf_material_archive()
        if initializer is None:
            pass
        elif type(initializer) == bytes:
            self._material_archive.ParseFromString(initializer)
        else:
            raise TypeError('Expected bytes. Found: ' + str(type(initializer)))

    def __len__(self) -> int:
        """Return the number of Material instances stored in the collection."""
        return len(self._material_archive.serialized_materials)

    def __iter__(self) -> MaterialArchiveIterator:
        """Iterate over Material objects stored in this collection.

        Preserves order of insertion.
        """
        return MaterialArchiveIterator(self)

    def __getitem__(self, index: int) -> Material:
        """Return the nth Material object in the collection.

        Preserves order of insertion.
        """
        try:
            index = int(index)
        except:
            raise TypeError(
                'Expected type that could be cast to int. Found: ' + str(type(index)))
        return Material(self._material_archive.serialized_materials[index])

    def append(self, material: Material) -> None:
        """Copies a Material object into the end of the collection.

        Increments length of the collection by one. Changes made to an
        object after it is passed to this method will not be reflected
        in the copy stored in the MaterialArchive.
        """
        if not isinstance(material, Material):
            raise TypeError(
                'Expected instance of Material. Found: ' + str(type(material)))
        self._material_archive.serialized_materials.append(
            material.serialize())
        
    def concatenate(self, other: MaterialArchive) -> MaterialArchive:
        '''
        Concatenate this MaterialArchive with another, returning the result.

        Repeated or shared elements will only appear in the output once.
        The order of the elements in the output may not be preserved.
        Neither of the input MaterialArchive objects are mutated.
        '''
        if not isinstance(other, MaterialArchive):
            raise TypeError("Expected MaterialArchive instance. Found: " + str(type(other)))
        bytes_set = set()
        for material_bytes in self._material_archive:
            bytes_set.add(material_bytes)
        for material_bytes in other._material_archive:
            bytes_set.add(material_bytes)
        output = MaterialArchive()
        for material_bytes in bytes_set:
            output._material_archive.append(output)
        return output

    def serialize(self) -> bytes:
        """Return a serialized version the wrapped protocol-buffer."""
        return self._material_archive.SerializeToString()

    def save_to_file(self, filename: str) -> None:
        """Serialize the wrapped protocol-buffer and save the result to a
        file."""
        if type(filename) != str:
            raise TypeError('Expected str. Found: ' + str(type(filename)))
        else:
            with open(filename, 'wb') as file:
                file.write(self.serialize())

    @staticmethod
    def load_from_file(filename) -> MaterialArchive:
        """Read a file, interpret the contents as a serialized MaterialArchive
        protocol-buffer, and return an instance of MaterialArchive wrapping the
        de- serialized protocol-buffer."""
        if type(filename) != str:
            raise TypeError('Expected str. Found: ' + str(type(filename)))
        else:
            with open(filename, 'rb') as file:
                return MaterialArchive(file.read())

    def to_numpy_double_array(self) -> numpy.ndarray:
        """Create a 2D numpy array from the application of
        Material.to_numpy_double_array to each of the Material objects stored
        in the collection.

        Axis 0 will have the same size as the length of the
        MaterialArchive collection. Axis 1 will have the same size as
        the length of the return value of
        Material.to_numpy_double_array.
        """
        rows = []
        for material in self:
            rows.append(material.to_numpy_double_array())
        return numpy.array(rows)

    def save_as_numpy(self, filename: str) -> None:
        """Convert the collection to a 2D numpy array using
        MaterialArchive.to_numpy_double_array and save the result in a npy
        file."""
        if type(filename) != str:
            raise TypeError('Expected str. Found: ' + str(type(filename)))
        else:
            numpy_double_array = self.to_numpy_double_array()
            numpy.save(filename, numpy_double_array,
                       allow_pickle=False, fix_imports=False)
