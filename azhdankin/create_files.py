#!python
""" File creation utility.

This is a utility to crate the directory and sample files which
will be used for transfer to a cloud storage.

The utility takes one command line parameter: number of files to create.
If the parameter is not given by default it will create 10 files.

The file creation is performed by copying the content of the ./seed-file.txt
content n-times (where n is a randomly generated number in a range from 1 to a 100)
into a destination file and naming the destination file by appending sequentially
incremented number to the base file name.

"""


import sys
import os
import random

#Set the root of the files location and the name prefix for the files to be generated
path = "./files/"
name_prefix = "file-2-upload"

#Set the default number of files to be generated
num_files = 10

#Read the number of files to be generated from the cmd line if provided
if len(sys.argv) > 1:
 num_files = int(sys.argv[1])

#Specify the "seed" for the generated files' content.
seed_file = "./seed-file.txt"

#Create the destination directory if it does not exist
if not os.path.exists(path):
  os.makedirs(path)

#Populate the seed string for the files to be created and initialize the content
file=open(seed_file,"r")
seed_content = file.read()
target_file_content = ""

#Create the files for upload
for target_file_idx  in range (0, num_files):
  #Replicate the seed content a random number of times
  repeat = random.randint(1,100)
  for chunk_num in range (0, repeat):
    target_file_content = target_file_content + seed_content 
  target_file = open (path + name_prefix + str(target_file_idx) + ".txt", 'w')
  target_file.write (target_file_content)  
  target_file_content = ""


