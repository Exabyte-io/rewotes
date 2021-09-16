import gevent
from gevent import monkey  # pylint: disable=import-error
from gevent.pool import Pool  # pylint: disable=import-error
monkey.patch_all() # patches blocking http I/O by non-blocking variant

import json, os, math
import time
import datetime
import libcloud
from libcloud.storage.types import Provider
from libcloud.storage.providers import get_driver
# print(libcloud.__file__)

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
        S = current_milli_time()
        # print('S', S)
        out = func(*args)
        E = current_milli_time()
        # print('E', E)
        print(str(func.__name__)+':', E-S, 'ms')
        return out
    return run

class BaseUpload:
    def __init__(self, provider, key, secret, container_name):
        self.provider = provider
        self.key = key
        self.secret = secret
        self.container_name = container_name
        self.factory = get_driver(provider)
        self.driver = self.initDriver()
        self.container = self.initContainer(container_name)

    def initDriver(self):
        return self.factory(key=self.key, secret=self.secret, signature_version='4')
    
    def uploadSimple(self, driver, src, dst, meta_data={}):
        extra = {'meta_data': meta_data} # {'acl': 'public-read'}
        obj = driver.upload_object(src, self.container, dst, extra=extra)
        return obj

    def uploadMultiPart(self, driver, src, dst, meta_data={}):
        extra = {'meta_data': meta_data} # {'acl': 'public-read'}
        with open(src, 'rb') as iterator:
            obj = driver.upload_object_via_stream(
                iterator=iterator,
                container=self.container,
                object_name=dst,
                extra=extra
            )
            return obj

    def initContainer(self, name):
        containers = self.driver.list_containers()
        container = next((el for el in containers if el.name==name), None)
        if not container:
            print('creating container...')
            container = self.driver.create_container(name)
        return container

    def upload(self, src, dst, meta_data={}, filesize=0):
        if (filesize/1024/1024 >= 50):
            return self.uploadMultiPart(self.driver, src, dst, meta_data)
        else:
            return self.uploadSimple(self.driver, src, dst, meta_data)

    def saveMeta(self):
        it = self.driver.iterate_container_objects(self.container)
        meta = {}
        for x in it:
            meta[x.name] = {
                'size': x.size,
                'modified': x.extra['last_modified'],
                'acl': None,
            }

        with open('meta.json', 'w') as outfile:
            json.dump(meta, outfile, indent=4)


class ParallelUpload(BaseUpload):
    def makeDummyFile(self, size_mb):
        f = open('dummy.dat',"wb")
        f.seek(1024*1024*size_mb-1)
        f.write(b"a")
        f.close()

    # @timed
    def asyncUpload(self, file, filesize=0):
        try:
            driver = self.initDriver() # new driver for each thread because not thread-safe
            meta = {'dt': str(datetime.datetime.now(datetime.timezone.utc)) }
            if (filesize/1024/1024 >= 50):
                obj = self.uploadMultiPart(driver, file, file, meta)
            else:
                obj = self.uploadSimple(driver, file, file, meta)
            return {'success': obj.name}
        except Exception as ex:
            # print(ex)
            return {'error': file}

    def findOptimalNoThreads(self, files):
        print('findOptimalNoThreads ...')
        # strategy: estimate upload speed, then divide by number of files to be uploaded
        # max 200 concurrent uploads
        hard_thread_limit = 200 # prevent too many threads

        avg_file_size = sum([f['size'] for f in files])/len(files)
        dummySizeMb = 5
        self.makeDummyFile(dummySizeMb)
        S = current_milli_time() # start time upload
        self.asyncUpload('dummy.dat', 10*1024*1024)
        E = current_milli_time() # end time upload
        dummyUploadSpeed = (E-S)/1000 # upload time in seconds
        upSpeed = dummySizeMb*1024*1024 / dummyUploadSpeed # bytes per second
        arbit_margin = 0.9 # arbitrary margin for misc. networking
        optimalThreads = max(1, math.floor(upSpeed / avg_file_size * arbit_margin)) # bytes/second / avg_bytes
        optimalThreads = min(hard_thread_limit, optimalThreads) 
        optimalThreads = min(len(files), optimalThreads) # don't create more threads than no. files
        return optimalThreads

    @timed
    def uploadFiles(self, files, poolsize):
        print('poolsize:', poolsize)
        pool = Pool(poolsize)
        tasks = []

        def spawn(tasks, path, filesize):
            g = pool.spawn(self.asyncUpload, path, filesize)
            g.meta = path;
            tasks += [g]

        for file in files:
            spawn(tasks, file['path'], file['size'])

        for g in tasks:
            try:
                out = g.get(timeout=5) # large files need large timeout (or none)
                if 'error' in out:
                    print('a.respawning', g.meta) # re-try mechanism
                    spawn(tasks, g.meta)
            except gevent.timeout.Timeout as ex:
                # print(ex)
                print('b.respawning', g.meta) # re-try mechanism
                spawn(tasks, g.meta)

