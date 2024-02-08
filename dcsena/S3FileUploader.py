from botocore.client import BaseClient

from dcsena.BaseFileUploader import BaseFileUploader


class S3FileUploader(BaseFileUploader):
    """S3 Implementation for FileUploader

    Args:
        bucket_name: S3 Bucket to upload files to
        s3_client: S3 Client
    """
    def __init__(self, bucket_name: str, s3_client: BaseClient):
        self.bucket_name = bucket_name
        self.s3_client = s3_client

    def upload_file(self, key: str, file_path: str):
        """Method for file upload

        Args:
            key: Key name of file to upload
            file_path: Path to the file to upload
        """
        # This function call is the same as S3Transfer upload_file
        # and is under the hood automatically doing a multipart upload when file is over a specific size.
        # This means that we don't need to handle specific logic for large files that would need to be chunked.
        # TransferConfig can be passed in if config changes are needed, such as concurrency and size limits
        self.s3_client.upload_file(file_path, self.bucket_name, key)
