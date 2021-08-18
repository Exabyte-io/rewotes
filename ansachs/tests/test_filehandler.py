import unittest

from filehandler import list_files, FileDesc


class FileHandlerTests(unittest.TestCase):

    def setUp(self) -> None:
        self.test_dir = 'testdir/'


    def test_list_files(self):
        files = list_files(self.test_dir)
        self.assertEqual(len(files), 3)
        assert FileDesc(path='testdir/nesteddir/file3.txt', filename='file3.txt') in files


if __name__ == '__main__':
    unittest.main()