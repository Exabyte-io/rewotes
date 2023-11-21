import sys
import os
import random

path = "./files/"
name_prefix = "file-2-upload"

num_files = 10

if len(sys.argv) > 1:
 num_files = int(sys.argv[1])

seed_file = "./seed-file.txt"

file=open(seed_file,"r")
seed_content = file.read()
target_file_content = ""

for target_file_idx  in range (0, num_files):
  repeat = random.randint(1,100)
  for chunk_num in range (0, repeat):
    target_file_content = target_file_content + seed_content 
  target_file = open (path + name_prefix + str(target_file_idx) + ".txt", 'w')
  target_file.write (target_file_content)  
  target_file_content = ""


