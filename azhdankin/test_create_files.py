#Test for file creation utility.
import glob
import os
from create_files import create_files

#Perform test of file creation
def test_create_files():

  #Set the root of the files location and the name prefix for the files to be generated
  path = "./test_files/"

  #Set the default number of files to be generated
  num_files = 10

  create_files(path, num_files)

  #Path to the root directory where the created files  are located
  path = "./test_files/*"

  #Initialize the file names list
  file_list = []

  #Populate the list of the file names that were generated. Currently we support two levels of the directories
  #where files are located
  for entry in glob.iglob(path, recursive=True):
    if os.path.isfile(entry):
      file_list.append(entry)
    else:
      entry = entry + "/*"
      for element in glob.iglob(entry, recursive=True):
        if os.path.isfile(element):
          file_list.append(element)

  file_list_len = len(file_list)
  assert file_list_len == 10
 




