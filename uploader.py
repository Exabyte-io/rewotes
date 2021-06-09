import boto3
from botocore.exceptions import ClientError
import logging
import threading

class Uploader:
    def __init__(self, client):
        self.client = client

    def upload_batch(self, buffer):
        jobs = []
        for file in buffer.files():
            target = lambda: self.upload_single_file(file)
            thread = threading.Thread(target=target)
            jobs.append(thread)
            thread.start()

        # Wait for all jobs to complete.
        for thread in jobs:
            thread.join()
        self.client.save_metadata()

    def upload_single_file(self, file):
        """ Uploads a single file through the client.
        """
        try:
            response = self.client.upload(file)
        except ClientError as e:
            logging.error(e)
