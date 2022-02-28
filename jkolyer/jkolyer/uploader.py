import os
import boto3
from botocore.exceptions import ClientError
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

class Uploader:
    def __init__(self):
        pass

    def get_uploaded_file(self, bucket_name, fname):
        client = boto3.client("s3")
        response = client.get_object(Bucket=bucket_name, Key=fname)
        contents = response["Body"].read()
        return contents

    def perform_upload(self):
        s3 = boto3.client('s3', region_name='us-east-1')
        s3.put_object(Bucket='mybucket', Key=self.name, Body=self.value)
                
    def upload_file(self, file_name, bucket, object_name=None):
        """Upload a file to an S3 bucket
        :param file_name: File to upload
        :param bucket: Bucket to upload to
        :param object_name: S3 object name. If not specified then file_name is used
        :return: True if file was uploaded, else False
        """
        # If S3 object_name was not specified, use file_name
        if object_name is None:
            object_name = os.path.basename(file_name)

        # Upload the file
        s3_client = boto3.client('s3')
        try:
            s3_client.upload_file(file_name, bucket, object_name)
        except ClientError as err:
            logging.error(err)
            return False
        
        return True

    
