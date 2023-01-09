import os
from mp_api.client import MPRester
from ..material import Material, MaterialArchive
from typing import Optional

class Downloader(object):
    """Manages downloads of materials from materialsproject.org."""

    def __init__(self, api_key: Optional[str]):
        """
        Create a Downloader instance.
        
        Argument is Materials Project API key encoded as a string, or None.
        If None, Downloader will load the API key from the "MP_API_KEY" environment variable.
        """
        if api_key is None:
            api_key = os.environ.get("MP_API_KEY")
        elif type(api_key) != str:
            raise TypeError('Expected str. Found: ' + str(type(api_key)))
        self.api_key = api_key

    def download(self, selected_fields: Optional[list] = None) -> MaterialArchive:
        """Download a selected list of material properties for all materials in
        the materialsproject.org database."""
        possible_field_list = [
            'composition',
            'composition_reduced'
        ]
        fields = [
            "material_id",
            'band_gap',
        ]
        if selected_fields is None:
            fields += possible_field_list
        elif type(selected_fields) == list:
            if len(selected_fields) == 0:
                raise ValueError("Expected nonempty list.")
            for field in selected_fields:
                if type(field) != str:
                    raise TypeError("Expected str. Found: " + str(type(field)))
                if field not in possible_field_list:
                    raise ValueError("Received invalid field name: " + field)
                if field in fields:
                    raise ValueError("Field already included in download list: " + field)
                fields.append(field)
        else:
            raise TypeError("Expected None or list. Found: " + str(type(selected_fields)))
        with MPRester(self.api_key) as mpr:
            material_ids = list(map(
                lambda document: document.material_id,
                mpr.summary.search(fields=['material_id'])
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
                        if "composition" in fields:
                            material.composition[atomic_number] = int(
                                document.composition[atomic_number])
                        if "composition_reduced" in fields:
                            material.composition_reduced[atomic_number] = int(
                                document.composition_reduced[atomic_number])
                    archive.append(material)
            return MaterialArchive

    def download_to_file(self, filename: str, selected_fields: Optional[list] = None) -> None:
        """Download a selected list of material properties for all materials in
        the materialsproject.org database.

        Save the results in a provided file.
        """
        # verify that we can write to the file before we download
        if type(filename) != str:
            raise TypeError('Expected str. Found: ' + str(type(filename)))
        with open(filename, 'w') as file:
            file.write('')
        archive = self.download(selected_fields)
        archive.save_to_file(filename)
