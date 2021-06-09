import os
import hashlib

BUF_SIZE = 65536

class File:
    def __init__(self, name):
        self.obj = open(name, 'rb')
        self.name = name
        self.size = os.path.getsize(name)
        self._hash = None

    def hash(self):
        """ Returns the hash of the file.
        """
        if self._hash is not None:
            return self._hash
        md5 = hashlib.md5()

        # Read the file object in chunks of BUF_SIZE.
        self.obj.seek(0)
        data = self.next()
        while data:
            md5.update(data)
            data = self.next()

        # Reset the file object to the beginning for upload.
        self.obj.seek(0)
        self._hash = md5.hexdigest()
        return self._hash

    def next(self):
        return self.obj.read(BUF_SIZE)

