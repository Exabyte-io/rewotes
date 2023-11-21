from google.cloud import storage

class CloudStorage:
    def __init__(self, name):
        self.name = name 

    def upload_object(self, object_name, source_file_name):
        pass
    
    def __str__(self):
        return self.name

class CloudStorageGCP(CloudStorage):
    def __init__(self, bucket_name, project=None):
        super().__init__("CloudStorageGCP")
        self.project = project
        self.bucket_name = bucket_name
        self.client = storage.Client(project=self.project)
        self.bucket = self.client.bucket(self.bucket_name)
        if not self.bucket.exists():
          self.bucket = self.client.create_bucket(self.bucket_name)

    def upload_object(self, object_name, source_file_name):
        # Create a new blob object
        blob = self.bucket.blob(object_name)
        # Upload the file to the bucket
        blob.upload_from_filename(source_file_name)


