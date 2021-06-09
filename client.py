import json
import datetime

from file import File

class Client:
    def __init__(self, id=0, key=0, bucket="main", is_testing=False):
        self.bucket = bucket
        self.is_testing = is_testing
        self.user = id

        if not self.is_testing:
            self.login(id, key)
            self.s3_client = boto3.client('s3')
        self.timestamp = str(datetime.datetime.now())
        self._metadata = None


    def login(self, id, key):
        """ Creats the credential folder for AWS.
        """
        with open("aws.config", 'r') as f:
            config = f.read()
        config.format(id=id, key=key)
        with open("~/.aws/credentials", 'w+') as f:
            f.write(config)


    def metadata(self):
        """ Returns the file that contains all the metadata.
        """
        if self._metadata is None:
            if self.is_testing:
                self._metadata = dict()
            else:
                self.s3_client.download_file(self.bucket,
                                             '/.metadata',
                                             '.metadata')
                with open('.metadata', 'r') as f:
                    self._metadata = json.load(f)
        return self._metadata

    def save_metadata(self):
        """ Uploads the updated metadata to the S3 bucket.
        """
        with open('.metadata', 'w') as f:
            f.write(json.dumps(self._metadata))
        self.upload(File('.metadata'))


    def upload(self, file):
        """ Uploads file to the S3 bucket.
        """
        if self.is_testing:
            print("Uploading %s..." % file.name)
        else:
            self.s3_client.upload_fileobj(file.obj, self.bucket, file.name)