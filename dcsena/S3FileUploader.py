from botocore.client import BaseClient

from dcsena.BaseFileUploader import BaseFileUploader


# S3 Implementation
class S3FileUploader(BaseFileUploader):

    def __init__(self, bucket_name: str, s3_client: BaseClient):
        self.bucket_name = bucket_name
        self.s3_client = s3_client

    def upload_file(self, key, file_name):
        # Same as S3 Transfer upload_file
        # Ths is under the hood automatically doing a multipart upload when file is over a specific size
        # TransferConfig can be passed in if config changes are needed, such as concurrency and size limits
        self.s3_client.upload_file(file_name, self.bucket_name, key)
