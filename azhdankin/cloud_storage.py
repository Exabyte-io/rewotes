"""
This module contains the definitions for the CloudStorage base class which acta like an "abstract class" or
the equivalent of the "interface" for the Cloud storage.
It also has an implementation of the Cloud storage implememntation for GCP, which allows to upload the files.
"""

from google.cloud import storage

#Base CloudStorage class
class CloudStorage:
    def __init__(self, name):
        self.name = name 

    def upload_object(self, object_name, source_file_name):
        pass
    
    def __str__(self):
        return self.name

#Cloud storage implementation for GCP
class CloudStorageGCP(CloudStorage):
    def __init__(self, bucket_name, project=None):
        super().__init__("CloudStorageGCP")
        self.project = project
        self.bucket_name = bucket_name
        self.client = storage.Client(project=self.project)

        #Resolve the reference to the destination bucket
        self.bucket = self.client.bucket(self.bucket_name)

        #If target bucket does not exist it will be created
        if not self.bucket.exists():
          self.bucket = self.client.create_bucket(self.bucket_name)

    def upload_object(self, object_name, source_file_name):
        # Create a new blob object
        blob = self.bucket.blob(object_name)
        # Upload the file to the bucket
        blob.upload_from_filename(source_file_name)


