import os
from mp_api.client import MPRester
from ..material import Material, MaterialArchive
from typing import Optional
from .interfaces import DownloaderInterface


class Downloader(DownloaderInterface):
    """Manages downloads of materials from materialsproject.org."""

    def __init__(self, api_key: Optional[str]):
        """Create a Downloader instance.

        Argument is Materials Project API key encoded as a string, or
        None. If None, Downloader will load the API key from the
        "MP_API_KEY" environment variable.
        """
        if api_key is None:
            api_key = os.environ.get('MP_API_KEY')
        elif type(api_key) != str:
            raise TypeError('Expected str. Found: ' + str(type(api_key)))
        self.api_key = api_key
        self.selected_fields = None
        self.query = None

    def set_feature_list(self, features: list) -> None:
        '''
        Set the list of features to download.

        Parameter is a list of strings, where each string is the name of a feature.
        If not used, all supported features will be downloaded.
        '''
        if type(features) != list:
            raise TypeError("Expected list. Found: " + str(type(features)))
        if len(features) == 0:
            raise ValueError("Must specify at least one feature to download.")
        for feature in features:
            if type(feature) != str:
                raise TypeError("Expected str. Found: " + str(type(feature)))
            if len(feature) == 0:
                raise ValueError("Feature name cannot be empty string.")
        self.selected_fields = features
    
    def set_query_elements_at_least(self, elements: list) -> None:
        '''
        Only download materials containing at least the specified elements.

        Parameter is a nonempty list of elemental symbols encoded as str.
        '''
        if type(elements) != list:
            raise TypeError("Expected list. Found: " + str(type(elements)))
        for element in elements:
            if type(element) != str:
                raise TypeError("Expected str. Found: " + str(type(element)))
            if len(element) < 1 or len(element) > 2:
                raise ValueError("Element symbol must be zero or two characters.")
        self.query = ("elements", elements)
    
    def set_query_elements_exact(self, elements: list) -> None:
        '''
        Only download materials containing exactly the specified elements.

        Parameter is a nonempty list of elemental symbols encoded as str.
        '''
        if type(elements) != list:
            raise TypeError("Expected list. Found: " + str(type(elements)))
        for element in elements:
            if type(element) != str:
                raise TypeError("Expected str. Found: " + str(type(element)))
            if len(element) < 1 or len(element) > 2:
                raise ValueError("Element symbol must be zero or two characters.")
        self.query = ("chemsys", elements)

    def download(self) -> MaterialArchive:
        """Download a selected list of material properties for selected materials in
        the materialsproject.org database."""
        possible_field_list = [
            'composition',
            'composition_reduced',
            "density",
            "density_atomic",
            "structure"
        ]
        fields = [
            'material_id',
            'band_gap',
        ]
        if self.selected_fields is None:
            fields += possible_field_list
        else:
            for field in self.selected_fields:
                if type(field) != str:
                    raise TypeError('Expected str. Found: ' + str(type(field)))
                if field not in possible_field_list:
                    raise ValueError('Received invalid field name: ' + field)
                if field in fields:
                    raise ValueError(
                        'Field already included in download list: ' + field)
                fields.append(field)
        with MPRester(self.api_key) as mpr:
            if self.query is None:
                material_ids = mpr.summary.search(fields=['material_id'])
            elif self.query[0] == "elements":
                material_ids = mpr.summary.search(elements=self.query[1], fields=['material_id'])
            elif self.query[0] == "chemsys":
                material_ids = mpr.summary.search(chemsys="".join(self.query[1]), fields=['material_id'])
            else:
                raise Exception("Interal Error")
            material_ids = list(map(
                lambda document: document.material_id,
                material_ids
            ))
            archive = MaterialArchive()
            for material_ids_slice_index in range(0, len(material_ids), 1000):
                material_ids_slice = material_ids[material_ids_slice_index:material_ids_slice_index+1000]
                docs = mpr.summary.search(
                    material_ids=material_ids_slice,
                    fields=fields
                )
                for document in docs:
                    material = Material()
                    material.material_id = document.material_id
                    material.band_gap = document.band_gap
                    for atomic_number in range(1, 119):
                        if 'composition' in fields:
                            material.composition[atomic_number] = int(
                                document.composition[atomic_number])
                        if 'composition_reduced' in fields:
                            material.composition_reduced[atomic_number] = int(
                                document.composition_reduced[atomic_number])
                    if "density" in fields:
                        material.density = document.density
                    if "density_atomic" in fields:
                        material.density_atomic = document.density_atomic
                    if "structure" in fields:
                        material.lattice_a = document.structure.lattice.a
                        material.lattice_b = document.structure.lattice.b
                        material.lattice_c = document.structure.lattice.c
                        material.lattice_alpha = document.structure.lattice.alpha
                        material.lattice_beta = document.structure.lattice.beta
                        material.lattice_gamma = document.structure.lattice.gamma
                        material.lattice_volume = document.structure.lattice.volume
                    archive.append(material)
            return archive

    def download_to_file(self, filename: str) -> None:
        """Download a selected list of material properties for all materials in
        the materialsproject.org database.

        Save the results in a provided file.
        """
        # verify that we can write to the file before we download
        if type(filename) != str:
            raise TypeError('Expected str. Found: ' + str(type(filename)))
        with open(filename, 'w') as file:
            file.write('')
        archive = self.download()
        archive.save_to_file(filename)
