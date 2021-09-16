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
### Linear upload example
```py
from uploader import *

"""
  provider: of type libcloud.storage.types.Provider
  key: cloud storage API key
  secret: cloud storage API secret
  container_name: bucket name, will be created if does not exist
"""
@timed
def test_linear_upload(provider, key, secret, container_name):
    print('test_linear_upload')
    uploader = BaseUpload(provider, key, secret, container_name)
    for file in getFiles('data/')[:20]: # provide location & get first 20 files
        path = file['path']
        meta = {'dt': str(datetime.datetime.now(datetime.timezone.utc)) }
        obj = uploader.upload(path, path, meta, file['size'])
        # print(path)

    uploader.saveMeta()

if __name__ == "__main__":
  test_linear_upload(Provider.S3, 'your_key', 'your_secret', 'your_bucket_name')
```

**Output:**
```
test_linear_upload
test_linear_upload: 3420 ms
```

### Parallel upload example
```py
from uploader import *

"""
  provider: of type libcloud.storage.types.Provider
  key: cloud storage API key
  secret: cloud storage API secret
  container_name: bucket name, will be created if does not exist
"""
@timed
def test_parallel_upload(provider, key, secret, container_name):
    print('test_parallel_upload')
    files = getFiles('data/')[:20] # provide location & get first 20 files
    uploader = ParallelUpload(provider, key, secret, container_name)
    poolsize = 20
    # poolsize = uploader.findOptimalNoThreads(files) # this op may take several seconds
    uploader.uploadFiles(files, poolsize)

    uploader.saveMeta()

if __name__ == "__main__":
  test_parallel_upload(Provider.S3, 'your_key', 'your_secret', 'your_bucket_name')
```
**Output:**
```
test_parallel_upload
poolsize: 20
uploadFiles: 561 ms
test_parallel_upload: 1149 ms
```

---

**NOTES**
- `findOptimalNoThreads(files)`: finds an optimal number of threads to use, with respect to the number of files to be uploaded, and the hard thread limit (default: 200)
- the uploader has a custom built-in re-try mechanism, which keeps running until upload succeeds. This can be prone to infinite loops in case uploading fails due unforeseen reasons.
- `ParallelUpload.uploadFiles` can have a thread timeout (seconds, default: None), in case it ever gets stuck. However large files may need a large timeout value -- it can be tricky to programmatically determine a good timeout value.
- libcloud's Driver is [not thread-safe](https://libcloud.readthedocs.io/en/stable/other/using-libcloud-in-multithreaded-and-async-environments.html), which is why `ParallelUpload.asyncUpload()` inits a new driver for each upload/thread.
- libcloud's Driver provides two upload functions: `upload_object` and `upload_object_via_stream`, the latter uses the multipart strategy whenever possible. However this is a slower method for small files. Which is why our code defaults to the regular upload strategy for files under 500MB.
