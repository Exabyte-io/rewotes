import unittest
import numpy
import torch
from pkg import DataRangeEncoder, DataRangeEncoderArray, Material, MaterialArchive, TrainingDataManager


class TestMaterial(unittest.TestCase):

    def test_material_can_be_instantiated(self):
        a = Material()

    def test_material_can_be_instantiated_from_empty_byte_string(self):
        a = Material(b'')

    def test_material_cannot_be_instantiated_from_other_types(self):
        with self.assertRaises(TypeError):
            a = Material(0)

    def test_material_can_be_reconstructed_from_serialized_form(self):
        a = Material()
        a.material_id = 'HelloWorld'
        b = Material(a.serialize())
        self.assertEqual(b.material_id, 'HelloWorld')

    def test_to_numpy_double_array_returns_numpy_double_array(self):
        a = Material()
        b = a.to_numpy_double_array()
        self.assertEqual(type(b), numpy.ndarray)

    def test_material_id_is_string(self):
        a = Material()
        self.assertEqual(type(a.material_id), str)

    def test_band_gap_is_flaot(self):
        a = Material()
        self.assertEqual(type(a.band_gap), float)

    def test_composition_is_dict_like(self):
        a = Material()
        a.composition[8] = 2
        self.assertEqual(a.composition[8], 2)

    def test_composition_reduced_is_dict_like(self):
        a = Material()
        a.composition_reduced[8] = 2
        self.assertEqual(a.composition_reduced[8], 2)


class TestMaterialArchive(unittest.TestCase):

    def test_material_archive_can_be_instantiated(self):
        a = MaterialArchive()

    def test_material_archive_can_be_instantiated_from_empty_byte_string(self):
        a = MaterialArchive(b'')

    def test_material_archive_cannot_be_instantiated_from_other_types(self):
        with self.assertRaises(TypeError):
            a = MaterialArchive(0)

    def test_material_archive_can_be_reconstructed_from_serialized_form(self):
        a = MaterialArchive()
        b = Material()
        b.material_id = 'B'
        a.append(b)
        c = MaterialArchive(a.serialize())
        self.assertEqual(len(c), 1)
        self.assertEqual(c[0].material_id, 'B')

    def test_to_numpy_double_array_returns_2d_numpy_double_array(self):
        a = MaterialArchive()
        b = Material()
        a.append(b)
        c = a.to_numpy_double_array()
        self.assertEqual(type(c), numpy.ndarray)
        self.assertEqual(len(c.shape), 2)

    def test_appending_material_increments_length(self):
        a = MaterialArchive()
        initial_length = len(a)
        a.append(Material())
        final_length = len(a)
        self.assertEqual(final_length, initial_length + 1)

    def test_can_iterate_over_material_archive(self):
        a = MaterialArchive()
        b = Material()
        b.material_id = 'B'
        a.append(b)
        loop_count = 0
        for material in a:
            loop_count += 1
            self.assertEqual(material.material_id, 'B')
        self.assertEqual(loop_count, 1)


class TestDataRangeEncoder(unittest.TestCase):

    def test_can_instantiate_data_range_encoder(self):
        a = DataRangeEncoder(5, 8)

    def test_encoded_value_is_between_zero_and_one(self):
        a = DataRangeEncoder(5, 8)
        b = a.encode(6)
        self.assertGreaterEqual(b, 0)
        self.assertLessEqual(b, 1)

    def test_encoded_value_roundtrips_to_same_value(self):
        a = DataRangeEncoder(5, 8)
        b = a.encode(6)
        c = a.decode(b)
        self.assertLess(abs(c-6), 10**-6)


class TestDataRangeEncoderArray(unittest.TestCase):

    def test_can_instantiate_data_range_encoder_array(self):
        a = DataRangeEncoderArray(numpy.array(
            [[1, 2], [3, 4], [5, 6]], dtype=numpy.double))

    def test_can_encode_array_of_same_length(self):
        a = DataRangeEncoderArray(numpy.array(
            [[1, 2], [3, 4], [5, 6]], dtype=numpy.double))
        b = a.encode(numpy.array([2, 5], dtype=numpy.double))
        self.assertEqual(type(b), numpy.ndarray)
        self.assertEqual(b.shape, (2,))
        self.assertEqual(b[0], 0.25)
        self.assertEqual(b[1], 0.75)

    def test_can_decode_array_of_same_length(self):
        a = DataRangeEncoderArray(numpy.array(
            [[1, 2], [3, 4], [5, 6]], dtype=numpy.double))
        b = a.decode(numpy.array([0, 1], dtype=numpy.double))
        self.assertEqual(type(b), numpy.ndarray)
        self.assertEqual(b.shape, (2,))
        self.assertEqual(b[0], 1)
        self.assertEqual(b[1], 6)


class TestTrainingDataManager(unittest.TestCase):

    def test_dataset_can_be_instantiated(self):
        a = MaterialArchive()
        a.append(Material())
        a.append(Material())
        a.append(Material())
        b = TrainingDataManager(a)
        self.assertEqual(type(b.training), torch.utils.data.DataLoader)
        self.assertEqual(type(b.testing), torch.utils.data.DataLoader)


if __name__ == '__main__':
    unittest.main()
