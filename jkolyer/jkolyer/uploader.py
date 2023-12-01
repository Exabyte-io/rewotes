"""Uploader: abstract superclass for wrapping object storage providers.
   
"""
import os
from abc import ABC, abstractmethod
import json
import boto3
from botocore.exceptions import ClientError
import logging
from moto import mock_s3  # workaround for multiprocessing / pytest limits

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)


class Uploader(ABC):
    """ static bucket_name """
    bucket_name = 'rewotes-pfu-bucket'

    @abstractmethod
    def get_uploaded_data(self, bucket_name, fname):
        """Retrieves stored data (either file or metadata) for given key.
           Used for testing and validation purposes.
        :param bucket_name: bucket where the data was uploaded
        :param key: lookup key for the uploaded data
        :return: bytes
        """
        pass

    @abstractmethod
    def upload_metadata(self, metadata, bucket_name, key):
        """Performs metadata upload to given bucket under given key.
        :param metadata: JSON string representation
        :param bucket_name: bucket where the data was uploaded
        :param key: lookup key for the uploaded data
        :return: bool: True if no errors, False otherwise
        """
        pass
                
    @abstractmethod
    def upload_file(self, file_name, bucket, key):
        """Upload a file to a given bucket
        :param file_name: File path to upload
        :param bucket: Bucket to upload to
        :param object_id: S3 object name
        :return: True if file was uploaded, else False
        """
        pass

    
class S3Uploader(Uploader):
    """S3Uploader: concrete instance for S3 object storage uploads.
       Uses `boto3` module.
    :Properties:  
    : client: instance provided by `boto3` for S3
    """
    """ process-based `boto3` client, which will vary based on environment and parallel algorithm"""
    boto3_client = None
    """ used to initialze `boto3` clients"""
    endpoint_url = None
    """ In multiprocessing tests, used to replace `endpoint_url` as a signal to use mock `boto3` client."""
    MOCK_ENDPOINT_KEY = 's3_mock'

    @classmethod
    def s3_mock(cls):
        """Initializes a mock wrapper for boto3 for use in testing purposes.
           To support tests with multiprocessing, we need to define the mock
           within this class (instead of keeping it with the tests).
           This may be candidate for refactoring.
        :return: mocked boto3 client
        """
        mock = mock_s3()
        mock.start()
        s3 = boto3.client('s3', region_name='us-east-1')
        s3.create_bucket(Bucket=cls.bucket_name)
        cls.endpoint_url = cls.MOCK_ENDPOINT_KEY
        return s3

    @classmethod
    def set_boto3_client(cls, client):
        """Sets the boto3 client used for uploads.  This is set one time, and 
           used by all instances of the receiver.  May be a mocked value in test.
        :param client: the `boto3` client used
        :return: None
        """
        cls.boto3_client = client
        client.create_bucket(Bucket=cls.bucket_name)
   
    @classmethod
    def set_endpoint_url(cls, url):
        cls.endpoint_url = url
   
    def __init__(self):
        """Instance constructor.  Sets `client` property.  
        """
        if self.boto3_client:
            ''' here the client was instantiated elsewhere '''
            self.client = self.boto3_client
        else:
            ''' we instantiate client here, based on endpoint_url '''
            if self.endpoint_url:
                ''' if we have endpoint that's not our mock value, use it in client '''
                if self.endpoint_url != self.MOCK_ENDPOINT_KEY:
                    self.client = boto3.client("s3", endpoint_url=self.endpoint_url)
                else:
                    ''' in this case we're running in test '''
                    self.client = self.s3_mock()
            else:
                ''' boto3 must be configured using environment variables '''
                self.client = boto3.client("s3")
        try:
            ''' make sure we have the expected bucket '''
            self.client.create_bucket(Bucket=self.bucket_name)
        except:
            logger.warn(f"Could not create bucket named {self.bucket_name}")
            self.client = None
    
    def get_uploaded_data(self, bucket_name, key):
        ''' see superclass '''
        response = self.client.get_object(Bucket=bucket_name, Key=key)
        contents = response["Body"].read()
        return contents

    def upload_metadata(self, metadata, bucket_name, key):
        ''' see superclass '''
        if not self.client:
            logger.warn(f"upload_metadadta:  no client for {key}")
            return False
        try:
            self.client.put_object(Bucket=bucket_name, Key=key, Body=json.dumps(metadata))
        except ClientError as err:
            logging.error(err)
            return False
        return True
                
    def upload_file(self, file_name, bucket, key, file_size):
        ''' see superclass '''
        if not self.client:
            logger.warn(f"upload_file:  no client for {key}: {file_name}")
            return False
        try:
            logger.info(f"S3Uploader.upload_file: {file_name}; file_size={file_size}")
            self.client.upload_file(file_name, bucket, key)
        except ClientError as err:
            logging.error(err)
            return False
        return True

    
