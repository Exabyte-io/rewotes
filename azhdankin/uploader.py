import os
from cloud_storage import CloudStorage

class FileUploader:
    def __init__(self, storage, files, start_idx, count):
        self.storage = storage 
        self.files = files
        self.start = start_idx
        self.count = count

    def run (self):
        for i in range(self.start, self.start + self.count):
          object_name = os.path.split(self.files[i])[1]
          self.storage.upload_object(object_name, self.files[i])


