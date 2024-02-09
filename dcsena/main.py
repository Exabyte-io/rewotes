#!/usr/bin/env python3
import concurrent.futures
import os
import boto3
import asyncio

from dcsena.AsyncDirectoryUploader import AsyncDirectoryUploader
from dcsena.S3FileUploader import S3FileUploader

BUCKET_NAME = os.environ["S3_BUCKET"]

RESOURCES_DIR_PATH = "./resources"


async def upload(dir_name):
    boto3_session = boto3.Session()
    s3_client = boto3_session.client('s3')
    s3_file_uploader = S3FileUploader(BUCKET_NAME, s3_client)
    executor = concurrent.futures.ThreadPoolExecutor()
    dir_uploader = AsyncDirectoryUploader(s3_file_uploader, executor)
    return await dir_uploader.upload_directory_async(dir_name)


if __name__ == "__main__":
    loop = asyncio.get_event_loop()
    loop.run_until_complete(upload(RESOURCES_DIR_PATH))
