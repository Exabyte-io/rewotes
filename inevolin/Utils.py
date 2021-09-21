import os
import time

def getFiles(dir):
    listOfFiles = list()
    for (path, dirnames, filenames) in os.walk(dir):
        listOfFiles += [{
                'path': os.path.join(path, file).replace('\\','/'),
                'name': file,
                'size': os.path.getsize(os.path.join(path, file))}
            for file in filenames]
    return listOfFiles

def current_milli_time():
    return round(time.time() * 1000)

def timed(func):
    def run(*args):
        start_ts = current_milli_time()
        out = func(*args)
        end_ts = current_milli_time()
        print(str(func.__name__)+':', end_ts-start_ts, 'ms')
        return out
    return run


def makeDummyFile(size_mb):
    f = open('dummy.dat',"wb")
    f.seek(1024*1024*size_mb-1)
    f.write(b"a")
    f.close()
