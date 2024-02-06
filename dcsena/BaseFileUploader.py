# from zope.interface import Interface
from abc import ABC


# Interface for File Uploading
class BaseFileUploader(ABC):
    def upload_file(self, key, file_name):
        pass
