# Multi-threaded file uploader

# Overview

This is a python application that uploads a set of given files to a cloud object storage in parallel through the cloud provider's or third party API.

# Features 

1. Support up to 100,000nds of files, all inside one directory with arbitrary sizes. The root directory may contain subdirectories.
2. The object storage container holds the objects is private and only credential-based access is allowed.
3. Each object inside object storage has an associated metadata which contains file size, last modification time and file permissions.

 The utility is fast (utilizes full network bandwidth), consumes low CPU (low enough not block all other processes) and low Memory (<25%)
 It supports GCP Cloud Storage, however it has a modular and Object oriented implementation so the other cloud providers can be added. 

# Prerequisites
  You must have Python 3.8 and the Google Cloud Storage Python client installed.

# To run
  1.  Clone the git repository. Make sure the create_files.py and upload_files.py file permissions are set to "executable".
  2.  Run ./create_files.py utility. This utility will create the files that need to be uploaded to the cloud storage. You can set the 
      number of file sto be created as a cmd line parameter (i.e. ./create_files.py 10000 to create 10000 files).
  3. Run ./upload_files.py to upload the files to the Cloud Storage.

# Notes
