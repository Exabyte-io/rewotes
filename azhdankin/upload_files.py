#!python
""" Main module to run to upload the files to the Cloud storage.

This program establishes the connection to the Cloud Storage (GCP in this case),
reads the names of the files available for upload
and performs parallel upload of the files to the specified Cloud Storage bucket.

"""
import os
import glob
import time

from concurrent.futures import ThreadPoolExecutor

from cloud_storage import CloudStorageGCP
from uploader import FileUploader

#Path to the root directory where the files to be uploaded are located
path = "./files/*"

#Initialize the file names list
file_list = []

#Populate the list of the file names set for upload. Currently we support two levels of the directories
#where files are located
for entry in glob.iglob(path, recursive=True):
  if os.path.isfile(entry):
    file_list.append(entry)
  else:
    entry = entry + "/*"
    for element in glob.iglob(entry, recursive=True):
      if os.path.isfile(element):
        file_list.append(element)

#Specify the maximum number of the workers that perform files upload simultaneously
MAX_UPLOAD_WORKERS = 100 

#Calculate the partitioning of the file names list - each partition or chunk will be assigned
#to a single upload worker

file_list_len = len(file_list)
step = int(file_list_len/MAX_UPLOAD_WORKERS)
remainder = file_list_len%MAX_UPLOAD_WORKERS

#Initialize a Cloud Storage Provider
storage = CloudStorageGCP("azhdanki-test-bucket1", project='rewotes')

#Create the Thread Pool which will be used to run the uploader tasks
pool = ThreadPoolExecutor (max_workers=MAX_UPLOAD_WORKERS)

#Schedule the upload tasks
i=0
time_start = time.time()

while i < (file_list_len - remainder):
  print(i)
  uploader = FileUploader (storage, file_list, i, step)
  pool.submit (uploader.run)
  i += step

if remainder > 0:
  uploader = FileUploader (storage, file_list, i, remainder)
  pool.submit (uploader.run)

pool.shutdown (wait=True)

time_end = time.time()
time_delta = time_end - time_start
print ("It took " + str(time_delta) + " seconds to upload " + str(file_list_len) + " files.")



