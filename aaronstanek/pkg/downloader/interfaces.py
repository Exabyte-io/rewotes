from abc import ABC
from ..material import MaterialArchiveInterface
from typing import Optional


class DownloaderInterface(ABC):
    def download(self) -> MaterialArchiveInterface:
        raise NotImplementedError()

    def download_to_file(self, filename: str) -> None:
        raise NotImplementedError()
