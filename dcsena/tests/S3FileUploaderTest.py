import unittest
from unittest.mock import Mock, patch

from dcsena.S3FileUploader import S3FileUploader

TEST_BUCKET = "foo123"


class S3FileUploaderTest(unittest.TestCase):
    @patch("botocore.client.BaseClient")
    def test_upload_happy(self, mock_s3_client):
        s3_file_uploader = S3FileUploader(TEST_BUCKET, mock_s3_client)
        s3_file_uploader.upload_file("my_key", "my_file")

    @patch("botocore.client.BaseClient")
    def test_upload_throws_error(self, mock_s3_client):
        mock_s3_client.upload_file.side_effect = Exception("Error")
        s3_file_uploader = S3FileUploader(TEST_BUCKET, mock_s3_client)
        with self.assertRaises(Exception):
            s3_file_uploader.upload_file("my_key", "my_file")


if __name__ == '__main__':
    unittest.main()
