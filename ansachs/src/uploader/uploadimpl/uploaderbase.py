from abc import ABC, abstractmethod


class UploaderBase(ABC):
    @abstractmethod
    def upload(self, file_path: str, file_key: str):
        pass
