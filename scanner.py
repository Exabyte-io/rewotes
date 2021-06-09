import os
import pyspeedtest

from buffer import Buffer
from uploader import Uploader

# The average portion of files that are updated per hour.
UPDATE_FREQUENCY = 0.1

class Scanner:

    def __init__(self, client, root_dir=None):
        if root_dir is None:
            root_dir = os.path.dirname(os.path.realpath(__file__))

        self.root_dir = root_dir
        self.client = client
        self.speedtest = pyspeedtest.SpeedTest("https://aws.amazon.com/")

    def sync(self):
        """ Checks each file in the root directory and uploads
        them through client.
        """
        bandwidth = self.speedtest.upload()
        chunk_size = bandwidth / UPDATE_FREQUENCY

        filenames = self.filenames()
        buffer = Buffer(self.client)
        uploader = Uploader(self.client)
        for file in filenames:
            buffer.add(file)

            # If we have reached our chunk size,
            # we clear out the buffer.
            if buffer.size() >= chunk_size:
                uploader.upload_batch(buffer)
                buffer.clear()

    def filenames(self):
        """ Returns a list of all file names in root directory
        and its sub-directories.
        """
        directory_list = os.listdir(self.root_dir)
        filenames = []

        while directory_list:
            entry = directory_list.pop(0)

            if os.path.isdir(entry):
                directory_list += [os.path.join(entry, name) \
                                        for name in os.listdir(entry)]
            else:
                filenames.append(entry)

        return filenames

