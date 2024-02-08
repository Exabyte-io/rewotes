from abc import ABC


class BaseFileUploader(ABC):
    """Abstract class for File Uploading"""
    def upload_file(self, key: str, file_path: str):
        """Abstract method for file upload

        Args:
            key: Key name of file to upload
            file_path: Path to the file to upload
        """
        pass
