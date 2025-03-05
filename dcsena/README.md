# Multi-threaded file uploader (Backend)

> Ideal candidate: skilled python developer with solid knowledge of cloud and distributed systems.

# Overview

Create a python application that uploads a set of given files to a cloud object storage in parallel through the cloud provider's or third party API.

# Requirements

1. Support up to 100,000nds of files, all inside one directory with arbitrary sizes. The root directory may contain subdirectories.
2. The object storage container which holds the objects is private and only credential-based access is allowed.
3. Each object inside object storage should have an associated metadata which contains file size, last modification time and file permissions.

# Expectations

- Fast (utilize full network bandwidth), low CPU (do not block all other processes) and low Memory (<25% tentatively) file uploader
- Support for AWS S3
- Modular and Object oriented implementation (to add other cloud providers)
- Clean and documented code
- Tests

# Dev Design

## Multithreading vs Multiprocessing vs Asyncio
Probably the most important decision for this implementation is how to handle the parallelization required for the parallel file uploader.
Given the constraints (low CPU, low memory, high IO operations, variable IO speed depending on file size), asyncio is the ideal choice here as it's meant for parallel, heavy IO ops.
However, boto3 and s3's upload_file are not async. There's a library called aiobotocore that attempts to create an async boto3 library, but it has pretty limited support (does not support S3Transfer, only put). It also has a hardpinned botocore dependencies, which will be a PITA from a dependency management size.
Luckily, aibotocore is not really needed, and we can leverage the asyncio running loop and run_in_executor with a thread pool to get around these limitations.



# Setup

Setup virtualenv and install dependencies
```commandline
python3 -m venv venv/
source venv/bin/activate
pip install -r requirements.txt
```

Setup AWS credentials if you have not already. BotoCredentials are best handled through environment variables.

Create resources file:
```commandline
python3 create_resources_script.py
```

Run file uploader:
```commandline
BUCKET_NAME="foo" python3 main.py
```
or if you don't have a default AWS profile setup:
```commandline
BUCKET_NAME="foo" AWS_ACCESS_KEY_ID="123" AWS_SECRET_ACCESS_KEY="abc" python3 main.py
```

# Things not done & Future Improvements
- I only tested with 1k small files here. Proper benchmarks to ensure desired memory and CPU profile at 100k would be a nice next step here.

## Requirements Check-In
1. By leveraging asyncio and threadpool, we can support upload of 100k files. Requirements met.
2. By leveraging a cloud file uploader such as S3, we are abstracting away authentication/authorization to cloud provider. However, S3, GCP, etc. all easily support public/private files with credential-based access. Requirements met.
3. By leveraging a cloud file uploader such as S3, we also get associated metadata for free. S3 will store file size, last modification time, and file perms.
With regards to last modification time and file permissions, I'm assuming that the cloud provider's metadata is sufficient here. If we wanted to store local modification time and local file perms, a metadata file would be needed. For example, if I upload a file at 12 PM, then at 1PM I modify it, and then at 2PM I upload it, the cloud provider would only say lastModificationTime is 2pm.
Local file permissions (RWX) is a can of worms, given different operating systems/users but if that was needed, similarly a metadata file would be needed. Marking requirements met under the listed assumptions.
