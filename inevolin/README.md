# Multi-threaded file uploader (Backend)

> Ideal candidate: skilled python developer with solid knowledge of cloud and distributed systems.

# Overview

Create a python application that uploads a set of given files to a cloud object storage in parallel through the cloud provider's or third party API.

# Requirements

1. Support up to 100,000nds of files, all inside one directory with arbitrary sizes. The root directory may contain subdirectories.
1. The object storage container which holds the objects is private and only credential-based access is allowed.
1. Each object inside object storage should have an associated metadata which contains file size, last modification time and file permissions.

# Expectations

- Fast (utilize full network bandwidth), low CPU (do not block all other processes) and low Memory (<25% tentatively) file uploader
- Support for AWS S3
- Modular and Object oriented implementation (to add other cloud providers)
- Clean and documented code
- Tests

# Timeline

We leave exact timing to the candidate. Must fit Within 5 days total.

# Notes

- we can provide temporary credentials to access AWS/Azure.

# Submission
The base code is located in `uploader.py`. The implementation facilitates two core features: a basic (linear) uploader and a parallel uploader. We use the `gevent` library for async networking, which also uses `greenlets` for concurrent programming (threading).

## Requirements
- Python 3.5+
- libcloud `pip install apache-libcloud`
- gevent `pip install gevent`

## Usage

## Config
To use this library you need to provide the API Access Key, Secret Key and Bucket/Container name. The simplest way is to make a copy of `Config.sample.py` and name it `Config.py`. In this Config file you can provide your credentials, but also override any default parameters. You can also include the config dict directly in your code instead, the default format and key/value pairs are as shown here:
```py
config = {   
    'KEY': 'change_me',
    'SECRET': 'change_me',
    'CONTAINER_NAME': 'change_me',

    'MULTIPART_SIZE_MB': 500, # minimum file size for using MultiPart technique
    'MAX_THREADS': 200, # max no. threads
    'DUMMY_SIZE_MB': 5, # size of dummy file for determining upload speed
    'DUMMY_ALGO_MARGIN': 0.9, # margin used for computing optimal number of threads [0.01, 1.0]
    'THREAD_TIMEOUT': None,
}
```

Now you can initialize one of the uploader classes as such:
```py
from libcloud.storage.types import Provider

# from Config import config
# or
# config = { ... }

provider = Provider.S3
uploader = ParallelUploader(provider, config)
```

### Linear upload example
```py
from BaseUploader import BaseUploader
from libcloud.storage.types import Provider
from Utils import timed, getFiles
import datetime

from Config import config

@timed
def test_linear_upload(provider):
    print('test_linear_upload')
    uploader = BaseUploader(provider, config)
    for file in getFiles('data/')[:20]: # provide location & get first 20 files
        path = file['path']
        meta = {'last_modified': str(datetime.datetime.now(datetime.timezone.utc)) }
        obj = uploader.upload(path, path, meta, file['size'])
    uploader.saveMeta()

if __name__ == "__main__":
  test_linear_upload(Provider.S3)
```

**Output:**
```
test_linear_upload
test_linear_upload: 3420 ms
```

### Parallel upload example
```py
from ParallelUploader import ParallelUploader
from libcloud.storage.types import Provider
from Utils import timed, getFiles
import datetime

from Config import config

@timed
def test_parallel_upload(provider):
    print('test_parallel_upload')
    files = getFiles('data/')[:20] # provide location & get first 20 files
    uploader = ParallelUploader(provider, config)
    # poolsize = 20
    poolsize = uploader.findOptimalNoThreads(files) # this op may take several seconds
    uploader.uploadFiles(files, poolsize)
    uploader.saveMeta()

if __name__ == "__main__":
  test_parallel_upload(Provider.S3)
```
**Output:**
```
test_parallel_upload
findOptimalNoThreads ...
poolsize: 20
uploadFiles: 846 ms
test_parallel_upload: 5256 ms
```

---

**NOTES**
- `findOptimalNoThreads(files)`: finds an optimal number of threads to use, with respect to the number of files to be uploaded, and the hard thread limit (default: 200)
- the uploader has a custom built-in re-try mechanism, which keeps running until upload succeeds. This can be prone to infinite loops in case uploading fails due unforeseen reasons. This can be improved by tracking & limiting retries per file.
- `ParallelUpload.uploadFiles` can have a thread timeout (seconds, default: None), in case it ever gets stuck or creeps. However large files may need a large timeout value -- it can be tricky to programmatically determine a good timeout value.
- libcloud's Driver is [not thread-safe](https://libcloud.readthedocs.io/en/stable/other/using-libcloud-in-multithreaded-and-async-environments.html), which is why `ParallelUpload.asyncUpload()` inits a new driver for each upload/thread.
- libcloud's Driver uses blocking I/O operations, to facilitate multithreading we use `gevent.monkey` to patch these with non-blocking variants, [read more](https://libcloud.readthedocs.io/en/stable/other/using-libcloud-in-multithreaded-and-async-environments.html).
- libcloud's Driver provides two upload functions: `upload_object` and `upload_object_via_stream`, the latter uses the multipart strategy whenever possible. However this is a slower method for small files. Which is why our code defaults to the regular upload strategy for files under 500MB.
