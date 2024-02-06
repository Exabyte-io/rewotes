import concurrent.futures
import glob
import asyncio
import functools

from dcsena import BaseFileUploader
from dcsena.BaseDirectoryUploader import BaseDirectoryUploader


class AsyncDirectoryUploader(BaseDirectoryUploader):
    def __init__(self, file_uploader: BaseFileUploader):
        self.file_uploader = file_uploader

    @staticmethod
    def get_files_to_upload(root_dir):
        return glob.glob("*", root_dir=root_dir, recursive=True)

    async def upload_directory_async(self, root_dir):
        loop = asyncio.get_running_loop()

        executor = concurrent.futures.ThreadPoolExecutor()
        tasks = []
        for f in self.get_files_to_upload(root_dir):
            key = f.strip(root_dir)
            print("RootDir: " + root_dir + " Key: " + key + " f: " + f)
            tasks.append(loop.run_in_executor(executor, functools.partial(self.file_uploader.upload_file, key=key, file_name=f)))
        completed, pending = await asyncio.wait(tasks)
        results = [t.result() for t in completed]
        return results
