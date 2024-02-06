import os
import boto3
import asyncio

from dcsena.AsyncDirectoryUploader import AsyncDirectoryUploader
from dcsena.S3FileUploader import S3FileUploader

BUCKET_NAME = os.environ["S3_BUCKET"]


async def upload():
    boto3_session = boto3.Session()
    s3_client = boto3_session.client('s3')
    s3_file_uploader = S3FileUploader(BUCKET_NAME, s3_client)
    dir_uploader = AsyncDirectoryUploader(s3_file_uploader)
    return await dir_uploader.upload_directory_async(".")


if __name__ == "__main__":
    loop = asyncio.get_event_loop()
    loop.run_until_complete(upload())
