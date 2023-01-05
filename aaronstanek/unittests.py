import unittest
import numpy
from pkg import Material, MaterialArchive

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
        a.material_id = "HelloWorld"
        b = Material(a.serialize())
        self.assertEqual(b.material_id, "HelloWorld")
    
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
        b.material_id = "B"
        a.append(b)
        c = MaterialArchive(a.serialize())
        self.assertEqual(len(c), 1)
        self.assertEqual(c[0].material_id, "B")
    
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
        b.material_id = "B"
        a.append(b)
        loop_count = 0
        for material in a:
            loop_count += 1
            self.assertEqual(material.material_id, "B")
        self.assertEqual(loop_count, 1)

if __name__ == '__main__':
    unittest.main()
