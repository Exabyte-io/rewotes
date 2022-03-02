import os
import pathlib
import argparse
import pathlib
import sys
import logging
import boto3

import asyncio
from jkolyer.models.batch_model import BatchJobModel, parallel_upload_files
from jkolyer.models.file_model import FileModel
from jkolyer.uploader import S3Uploader

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)
        
def parse_cmd_line_arguments():
    """Describe
    :param name: describe
    :param name: describe
    :return: type describe
    """
    parser = argparse.ArgumentParser(
        description="PFU: Parallel File Upload",
        epilog="Thanks for using the service!",
    )
    parser.add_argument(
        "--parallel",
        action='store_true',
        help="Runs the uploads in multiple processes (up to CPU count), default is concurrent.",
    )
    parser.add_argument(
        "--concurrent",
        action='store_true',
        help="Runs the uploads in a single process using asyncio (default).",
    )
    parser.add_argument(
        "--endpoint_url",
        nargs=1,
        action="store",
        metavar="ENDPOINT_URL",
        help="Endpoint for localstack S3 in the form http://localhost:4566",
    )
    parser.add_argument(
        "--root_dir",
        metavar="ROOT_DIR",
        action="store",
        required=True,
        help="Directory to load files for upload",
    )
    return parser.parse_args()

def perform_file_upload(parallel, root_dir):
    logger.info(f"initializing database, file root = {root_dir}")
    
    FileModel.create_tables()
    BatchJobModel.create_tables()
    
    batch = BatchJobModel.new_instance(root_dir)
    batch.generate_file_records()
    batch.reset_file_status()

    logger.info(f"performing upload")
    if parallel:
        parallel_upload_files(batch)
    else:
        loop = asyncio.get_event_loop()
        try:
            loop.run_until_complete(
                asyncio.gather(
                    batch.async_upload_files()
                ))
        finally:
            loop.close()

if __name__ == '__main__':
    args = parse_cmd_line_arguments()
    root_dir = pathlib.Path(args.root_dir)
    if not root_dir.is_dir():
        print("The specified root directory doesn't exist")
        sys.exit()
        
    concurrent = args.concurrent
    parallel = args.parallel
    if concurrent and parallel:
        parallel = False

    endpoint_url = args.endpoint_url
    if endpoint_url:
        client = boto3.client("s3", endpoint_url=endpoint_url[0], region_name='us-east-1')
        S3Uploader.set_boto3_client(client)
        S3Uploader.set_endpoint_url(endpoint_url[0])
        
    perform_file_upload(parallel, root_dir)
    

