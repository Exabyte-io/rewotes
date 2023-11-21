import os
import glob
from concurrent.futures import ThreadPoolExecutor

from cloud_storage import CloudStorageGCP
from uploader import FileUploader


path = "./files/*"

file_list = []

for entry in glob.iglob(path, recursive=True):
  if os.path.isfile(entry):
    file_list.append(entry)
  else:
    entry = entry + "/*"
    for element in glob.iglob(entry, recursive=True):
      if os.path.isfile(element):
        file_list.append(element)

MAX_UPLOAD_WORKERS = 2 

file_list_len = len(file_list)

step = int(file_list_len/MAX_UPLOAD_WORKERS)
remainder = file_list_len%MAX_UPLOAD_WORKERS

storage = CloudStorageGCP("azhdanki-test-bucket1", project='rewotes')

pool = ThreadPoolExecutor (max_workers=MAX_UPLOAD_WORKERS)

i=0

while i < (file_list_len - remainder):
  uploader = FileUploader (storage, file_list, i, step)
  pool.submit (uploader.run())
  i += step

if remainder > 0:
  uploader = FileUploader (storage, file_list, i, remainder)
  pool.submit (uploader.run())

pool.shutdown (wait=True)


