import unittest
from asyncio import Future
from concurrent.futures import Executor
from threading import Lock
from unittest.mock import patch
from unittest import IsolatedAsyncioTestCase

from dcsena.AsyncDirectoryUploader import AsyncDirectoryUploader

TEST_BUCKET = "foo123"


class AsyncDirectoryUploaderTest(unittest.TestCase):
    @patch("dcsena.BaseFileUploader")
    @patch("concurrent.futures.Executor")
    def test_get_files(self, mock_file_uploader, mock_executor):
        sut = AsyncDirectoryUploader(mock_file_uploader, mock_executor)
        files = sut.get_files_to_upload("../resources")
        # files exist and they don't include directories
        self.assertTrue(len(files) != 0)
        self.assertTrue("child-dir-0" not in files)


# DummyExecutor as it's easier to use than mocking
class DummyExecutor(Executor):

    def __init__(self):
        self._shutdown = False
        self._shutdownLock = Lock()

    def submit(self, fn, *args, **kwargs):
        with self._shutdownLock:
            if self._shutdown:
                raise RuntimeError('cannot schedule new futures after shutdown')

            f = Future()
            try:
                result = fn(*args, **kwargs)
            except BaseException as e:
                f.set_exception(e)
            else:
                f.set_result(result)

            return f

    def shutdown(self, wait=True, *, cancel_futures=False):
        with self._shutdownLock:
            self._shutdown = True


class AsyncDirectoryUploaderTestAsync(IsolatedAsyncioTestCase):
    @patch("dcsena.BaseFileUploader")
    async def test_upload_directory_async_happy(self, mock_file_uploader):
        mock_executor = DummyExecutor()
        sut = AsyncDirectoryUploader(mock_file_uploader, mock_executor)
        r = await sut.upload_directory_async("../resources")
        self.assertTrue(r)

    @patch("dcsena.BaseFileUploader")
    async def test_upload_directory_async_unhappy(self, mock_file_uploader):
        mock_file_uploader.upload_file.side_effect = Exception("Error")
        mock_executor = DummyExecutor()
        sut = AsyncDirectoryUploader(mock_file_uploader, mock_executor)
        with self.assertRaises(Exception):
            await sut.upload_directory_async("../resources")


if __name__ == '__main__':
    unittest.main()
