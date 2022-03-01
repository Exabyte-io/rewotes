"""Uploader: abstract superclass for wrapping object storage providers.
"""
import os
from abc import ABC, abstractmethod
import boto3
from botocore.exceptions import ClientError
import logging
from moto import mock_s3  # workaround for multiprocessing / pytest limits

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

class Uploader(ABC):
    
    def __init__(self):
        pass

    @abstractmethod
    def get_uploaded_data(self, bucket_name, fname):
        pass

    @abstractmethod
    def upload_metadata(self, metadata, bucket_name, name):
        pass
                
    @abstractmethod
    def upload_file(self, file_name, bucket, object_id):
        pass

    
class S3Uploader(Uploader):
    """S3Uploader: concrete instance for S3 object storage uploads.
       Uses `boto3` module.
    :Properties:  
    : client: instance provided by `boto3` for S3
    """

    def __init__(self, in_test=False):
        """Instance constructor.  Sets `client` property
        :param in_test: workaround for overriding default `boto3`
           in multiprocessing test scenarios
        """
        if in_test:
            s3_mock = mock_s3()
            s3_mock.start()
                
        self.client = boto3.client("s3")
    
    def get_uploaded_data(self, bucket_name, key):
        """Retrieves stored data (either file or metadata) for given key.
        :param bucket_name: bucket where the data was uploaded
        :param key: lookup key for the uploaded data
        :return: bytes
        """
        response = self.client.get_object(Bucket=bucket_name, Key=key)
        contents = response["Body"].read()
        return contents

    def upload_metadata(self, metadata, bucket_name, key):
        """Performs metadata upload to S3 bucket under given key.
        :param metadata: JSON string representation
        :param bucket_name: bucket where the data was uploaded
        :param key: lookup key for the uploaded data
        :return: bool: True if no errors, False otherwise
        """
        try:
            self.client.put_object(Bucket=bucket_name, Key=key, Body=metadata)
        except ClientError as err:
            logging.error(err)
            return False
        return True
                
    def upload_file(self, file_name, bucket, object_id):
        """Upload a file to an S3 bucket
        :param file_name: File path to upload
        :param bucket: Bucket to upload to
        :param object_id: S3 object name
        :return: True if file was uploaded, else False
        """
        try:
            self.client.upload_file(file_name, bucket, object_id)
        except ClientError as err:
            logging.error(err)
            return False
        
        return True

    
