"""
Class performing upload of the files to the cloud storage.
"""

import os

class FileUploader:
    """
    Constructor is taking the following parameters:
      storage - a reference to the instance of CloudStorage object.
                a CloudStorage class is a parent class for the cloud storage
                provider specific implementations, i.e. CloudStorageGCP, CloudStrageAWS

      files - a reference to the entire list of the file names that need to be uploaded

      start_idx - an index, a "pointer" to the files list specifying where this file uploader 
                  will start
 
      count - a number which specifies how many files this instance of uploader has to process (upload)

     """
    def __init__(self, storage, files, start_idx, count):
        self.storage = storage 
        self.files = files
        self.start = start_idx
        self.count = count

    #The method performing the upload of the group of the files to the Cloud storage
    def run (self):
        for i in range(self.start, self.start + self.count):
          object_name = os.path.split(self.files[i])[1]
          self.storage.upload_object(object_name, self.files[i])


