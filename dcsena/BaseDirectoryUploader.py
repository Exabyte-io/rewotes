from abc import ABC


# Interface for Directory Uploading
class BaseDirectoryUploader(ABC):
    def upload_directory_async(self, root_dir):
        pass
