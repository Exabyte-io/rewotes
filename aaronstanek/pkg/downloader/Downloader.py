import time
from mp_api.client import MPRester
from ..material import Material, MaterialArchive


class Downloader(object):
    def __init__(self, destination_filename: str, api_key: str):
        if type(destination_filename) != str:
            raise TypeError('Expected str. Found: ' +
                            str(type(destination_filename)))
        if type(api_key) != str:
            raise TypeError('Expected str. Found: ' + str(type(api_key)))
        if destination_filename[-16:] != '.materialarchive':
            raise ValueError(
                'Expected file extension materialarchive. Found filename: ' + destination_filename)
        # verify that we are able to write to the file before we begin downloading data
        with open(destination_filename, 'w') as file:
            file.write('')
        self.destination_filename = destination_filename
        self.api_key = api_key

    def download(self) -> MaterialArchive:
        start_time = time.time()
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
                    fields=[
                        'material_id',
                        'band_gap',
                        'composition',
                        'composition_reduced'
                    ]
                )
                for document in docs:
                    material = Material()
                    material.material_id = document.material_id
                    material.band_gap = document.band_gap
                    for atomic_number in range(1, 119):
                        material.composition[atomic_number] = int(
                            document.composition[atomic_number])
                        material.composition_reduced[atomic_number] = int(
                            document.composition_reduced[atomic_number])
                    archive.append(material)
            return MaterialArchive
