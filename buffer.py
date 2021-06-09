import gc
import logging

from file import File


class Buffer:
    """ Class that contains list of file objects and names.
    Responsible for updating metadata for each file.
    """
    def __init__(self, client):
        self.client = client
        self.metadata = self.client.metadata()
        self._size = 0
        self._files = []


    def add(self, filename):
        file = File(filename)
        self._size += file.size
        self._files.append(file)

    def clear(self):
        self._size = 0
        for file in self._files:
            del file.obj
        gc.collect()
        self._files = []

    def size(self):
        return self._size

    def files(self):
        files = []
        for file in self._files:
            if self.is_file_updated(file):
                self.update_metadata_for_file(file)
                files.append(file)

        return files

    def is_file_updated(self, file):
        """ Checks if file is different from that stored in the cloud.
        """
        user = self.client.user
        if file.name not in self.metadata:
            return True
        file_metadata = self.metadata.get(file.name, {})
        if user in file_metadata.get("permissions", (user, )):
            return file.hash() == file_metadata.get("hash", "")
        logging.error("Permission Denied for %s" % file.name)
        return False

    def update_metadata_for_file(self, file):
        """ Updates self.metadata for file.
        """
        file_metadata = dict()
        file_metadata["permissions"] = [self.client.user]
        file_metadata["hash"] = file.hash()
        file_metadata["last_updated"] = self.client.timestamp
        file_metadata["file_size"] = file.size
        self.metadata[file.name] = file_metadata