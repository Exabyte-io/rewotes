#!/usr/bin/env python3
import os
import boto3

BUCKET_NAME = os.environ["S3_BUCKET"]

RESOURCES_DIR_PATH = "./resources"


def validate_results(s3_client):
    paginator = s3_client.get_paginator('list_objects_v2')
    page_iterator = paginator.paginate(Bucket=BUCKET_NAME)
    num_files = 0
    for page in page_iterator:
        for file in page['Contents']:
            print(file['Key'])
            print(file['LastModified'])
            print(file['Size'])
            num_files += 1
    print(num_files)


if __name__ == "__main__":
    boto3_session = boto3.Session()
    s3_client = boto3_session.client('s3')
    validate_results(s3_client)
