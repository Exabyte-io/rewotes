from abc import ABC
from ..material import MaterialArchiveInterface
from typing import Optional


class DownloaderInterface(ABC):
    def download(self, selected_fields: Optional[list] = None) -> MaterialArchiveInterface:
        raise NotImplementedError()

    def download_to_file(self, filename: str, selected_fields: Optional[list] = None) -> None:
        raise NotImplementedError()
