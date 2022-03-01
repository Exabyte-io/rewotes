import os
from abc import ABC, abstractmethod
import boto3
from botocore.exceptions import ClientError
import logging

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

    def __init__(self):
        self.client = boto3.client("s3")
    
    def get_uploaded_data(self, bucket_name, fname):
        response = self.client.get_object(Bucket=bucket_name, Key=fname)
        contents = response["Body"].read()
        return contents

    def upload_metadata(self, metadata, bucket_name, name):
        try:
            self.client.put_object(Bucket=bucket_name, Key=name, Body=metadata)
        except ClientError as err:
            logging.error(err)
            return False
        return True
                
    def upload_file(self, file_name, bucket, object_id):
        """Upload a file to an S3 bucket
        :param file_name: File to upload
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

    
