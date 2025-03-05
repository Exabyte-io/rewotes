import os
import glob
import asyncio
import functools
from concurrent.futures import Executor

from dcsena import BaseFileUploader
from dcsena.BaseDirectoryUploader import BaseDirectoryUploader


class AsyncDirectoryUploader(BaseDirectoryUploader):
    """Class for uploading a directory using the fileuploader implementation for storage

    Args:
        file_uploader: Client for uploading files through whatever storage provider
        executor: thread pool executor to execute tasks in parallel
    """
    def __init__(self, file_uploader: BaseFileUploader, executor: Executor):
        self.file_uploader = file_uploader
        self.executor = executor

    @staticmethod
    def get_files_to_upload(root_dir: str):
        """Gets all files recursively in the specified directory

        Args:
            root_dir: root directory to begin file search

        Returns:
            All files in directory
        """
        files = glob.glob("**/*", root_dir=root_dir, recursive=True)
        # We can't upload directories so first need to filter them out
        return [f for f in files if os.path.isfile(os.path.join(root_dir, f))]

    async def upload_directory_async(self, root_dir: str):
        """Uses asyncio and executor to make async file upload calls for all files in directory

        Args:
            root_dir: root directory to upload files for

        Returns:
            completed task results
        """
        loop = asyncio.get_running_loop()
        tasks = []
        for f in self.get_files_to_upload(root_dir):
            file_path = os.path.join(root_dir, f)
            key = f.replace("\\", "/")
            tasks.append(loop.run_in_executor(self.executor, functools.partial(self.file_uploader.upload_file, key=key,
                                                                               file_path=file_path)))
        completed, pending = await asyncio.wait(tasks)
        results = [t.result() for t in completed]
        return results
