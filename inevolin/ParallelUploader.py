import math
import datetime
from gevent.pool import Pool  # pylint: disable=import-error
from BaseUploader import BaseUploader
from Utils import timed, current_milli_time, makeDummyFile

class ParallelUploader(BaseUploader):
    # @timed
    def asyncUpload(self, file, filesize=0):
        try:
            driver = self.initDriver() # new driver for each thread because not thread-safe
            meta = {'last_modified': str(datetime.datetime.now(datetime.timezone.utc)) }
            # use multipart upload for +X MB files
            if (filesize/1024/1024 >= self.config.get('MULTIPART_SIZE_MB')): 
                obj = self.uploadMultiPart(driver, file, file, meta)
            else:
                obj = self.uploadSimple(driver, file, file, meta)
            return {'success': obj.name}
        except Exception as ex:
            # print(ex)
            return {'error': file}

    def findOptimalNoThreads(self, files):
        # strategy: Estimate upload speed by uploading a dummy file
        #           then divide by number of files to be uploaded.
        print('findOptimalNoThreads ...')
        # max concurrent uploads, prevent too many threads
        hard_thread_limit = self.config.get('MAX_THREADS') #
        avg_file_size = sum([f['size'] for f in files]) / len(files)
        dummySizeMb = self.config.get('DUMMY_SIZE_MB')
        makeDummyFile(dummySizeMb) # makes a dummy file for uploading
        upload_time_start = current_milli_time() # start time upload
        self.asyncUpload('dummy.dat', 10*1024*1024)
        upload_time_end = current_milli_time() # end time upload
        dummyUploadSpeed = (upload_time_end-upload_time_start)/1000 # upload time in seconds
        upSpeed = dummySizeMb*1024*1024 / dummyUploadSpeed # bytes per second
        arbit_margin = self.config.get('DUMMY_ALGO_MARGIN') # arbitrary margin for misc. networking
        optimalThreads = max(1, math.floor(upSpeed / avg_file_size * arbit_margin)) # bytes/second / avg_bytes
        optimalThreads = min(hard_thread_limit, optimalThreads) 
        optimalThreads = min(len(files), optimalThreads) # don't create more threads than no. files
        return optimalThreads

    def spawnThread(self, pool, tasks, path, filesize):
        task = pool.spawn(self.asyncUpload, path, filesize)
        task.meta = path;
        tasks += [task]
    
    @timed
    def uploadFiles(self, files, poolsize):
        print('poolsize:', poolsize)
        pool = Pool(poolsize)
        tasks = []

        for file in files:
            self.spawnThread(pool, tasks, file['path'], file['size'])

        for task in tasks:
            try:
                out = task.get(timeout=self.config.get('THREAD_TIMEOUT')) # large files need large timeout (or none)
                if 'error' in out:
                    print('a.respawning', task.meta) # re-try mechanism
                    self.spawnThread(pool, tasks, task.meta)
            except gevent.timeout.Timeout as ex:
                # print(ex)
                print('b.respawning', task.meta) # re-try mechanism
                self.spawnThread(pool, tasks, task.meta)

