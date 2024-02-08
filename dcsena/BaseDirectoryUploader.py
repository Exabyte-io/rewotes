from abc import ABC


class BaseDirectoryUploader(ABC):
    def upload_directory_async(self, root_dir: str):
        """Abstract class for Directory Uploading

        Args:
            root_dir: Directory to upload
        """
        pass
