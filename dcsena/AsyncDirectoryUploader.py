import os
import glob
import asyncio
import functools
from concurrent.futures import Executor

from dcsena import BaseFileUploader
from dcsena.BaseDirectoryUploader import BaseDirectoryUploader


# Class for uploading a directory using the file uploader's implementation for storage
class AsyncDirectoryUploader(BaseDirectoryUploader):
    def __init__(self, file_uploader: BaseFileUploader, executor: Executor):
        self.file_uploader = file_uploader
        self.executor = executor

    @staticmethod
    def get_files_to_upload(root_dir: str):
        files = glob.glob("**/*", root_dir=root_dir, recursive=True)
        # We can't upload directories so first need to filter them out
        return [f for f in files if os.path.isfile(os.path.join(root_dir, f))]

    # Using asyncio and thread pool executor to make async file upload calls
    async def upload_directory_async(self, root_dir: str):
        loop = asyncio.get_running_loop()
        tasks = []
        for f in self.get_files_to_upload(root_dir):
            file_name = os.path.join(root_dir, f)
            key = f.replace("\\", "/")
            tasks.append(loop.run_in_executor(self.executor, functools.partial(self.file_uploader.upload_file, key=key,
                                                                               file_name=file_name)))
        completed, pending = await asyncio.wait(tasks)
        results = [t.result() for t in completed]
        return results
