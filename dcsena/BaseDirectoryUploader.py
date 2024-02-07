from abc import ABC


# Abstract class for Directory Uploading
class BaseDirectoryUploader(ABC):
    def upload_directory_async(self, root_dir: str):
        pass
