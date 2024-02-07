from abc import ABC


# Abstract class for File Uploading
class BaseFileUploader(ABC):
    def upload_file(self, key: str, file_name: str):
        pass
