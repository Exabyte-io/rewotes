# Parallel uploading service in python


## Objective

Create a python application that uploads a set of given files to a cloud object storage in parallel through the cloud provider's or third party API.

### Requirements

1. Support up to 100,000nds of files, all inside one directory with arbitrary sizes. The root directory may contain subdirectories.
1. The object storage container which holds the objects is private and only credential-based access is allowed.
1. Each object inside object storage should have an associated metadata which contains file size, last modification time and file permissions.

### Expectations

- Fast (utilize full network bandwidth), low CPU (do not block all other processes) and low Memory (<25% tentatively) file uploader
- Support for AWS S3
- Modular and Object oriented implementation (to add other cloud providers)
- Clean and documented code
- Tests

## Implementation

### Overview
This service ***pfu*** (parallel file upload) reads files from a root directory and loads a SQLite database with file metadata records.  Then a batch job initiates uploading files and metadata to a given object storage provider iteratively.

Parallism is achieved through [concurrency](https://docs.python.org/3/library/asyncio.html) or through [multiprocessing](https://docs.python.org/3/library/multiprocessing.html) – these are the two modes which the service runs.

The service has a test suite which mocks object storage.  In regular usage it leverages [localstack](https://localstack.cloud/) as storage provider.

Speed is achieved through parallelism, low CPU through using maximum `nice` setting, and low memory through paginating database result sets.

A companion file generator utility creates a tree of sample files, with user-specified depth.  This is for development and testing purposes only.

### Environment
* Python 3.9
* SQLite 3
* [mkfile](https://ss64.com/bash/mkfile.html) (sample tests)
* [localstack](https://localstack.cloud)
* Virtualenv
* Shell scripting

#### Inputs
* File root directory
* Localstack docker instance
* Object storage provider endpoint

#### Outputs
* Files and metadata uploaded into object storage
* SQLite database file
	* file metadata
	* upload result flag (success/fail)
	* batch upload metadata


## Install & Run

### Installation
### Run the service
Shell scripts are provided as a convenience.  The *file sample directory* is hardcoded as "samples".

* `cd $PROJECT_ROOT`
* `source scripts.sh`
* `activate` – for the virtualenv
* Create samples directory tree
	* `samplegen` script, providing tree depth argument; e.g., `samplegen 3` will create 3-level directory structure, each level containing 3 files (using `mkfile` utility).
* Run localstack docker instance
	* follow localstack instructions 
	* note which port it's running 
* Parallel mode
	* `localstack_p` with port number argument; e.g., `localstack_p 4566`
* Concurrent mode
	* `localstack_c` with port number argument; e.g., `localstack_p 4566`
* Tests use the `runtest` script.  

## Design Components

### Data model
#### Database Tables
##### FileStat
Tracks file metadata with upload status.  
> `CREATE TABLE FileStat
                  ( id TEXT PRIMARY KEY, 
                    created_at INTEGER,
                    file_size INTEGER,
                    last_modified INTEGER,
                    permissions TEXT,
                    file_path TEXT,
                    status INTEGER
                  )`

##### BatchJob
Tracks upload sessions and root directory.  
> `CREATE TABLE BatchJob
        ( id TEXT PRIMARY KEY, 
        status INTEGER,
        created_at INTEGER,
        root_dir TEXT
        )`

#### Data Objects
##### BaseModel
Abstract base class for other models.  Provides support for SQL generation.

##### FileModel
Represents file for upload, and it's upload status.  Provides SQL generation, and upload state machine.

##### BatchJobModel
Tracks and performs upload sessions.  Creates `FileModel` instances for a given directory and provides iteration across data.  Handles parallel and concurrent upload process.  SQL generation.  


### Uploader

The abstraction class `Uploader` serves as an interface for the `S3Uploader`.  It performs the basic unit of work for the service:  uploading a file and metadata to 3rd party provider.  Uploads are performed by order of ascending file size.

Upload jobs can be run in concurrent mode (using [asyncio](https://docs.python.org/3/library/asyncio.html)) or parallel mode (using [multiprocessing](https://docs.python.org/3/library/multiprocessing.html)).  

### File generator

Bootstrapping a test environment.  Creates a small directory tree with files of varying sizes (under 1G).  