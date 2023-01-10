from abc import ABC
from ..material import MaterialArchiveInterface
from typing import Optional


class DownloaderInterface(ABC):
    """Manages downloads of materials from a remote source."""
    def download(self) -> MaterialArchiveInterface:
        """Download a selected list of material properties for selected materials."""
        raise NotImplementedError()

    def download_to_file(self, filename: str) -> None:
        '''Download a selected list of material properties for selected materials, saving the result to a file.'''
        raise NotImplementedError()
