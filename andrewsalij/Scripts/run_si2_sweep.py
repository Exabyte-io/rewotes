import andrewsalij.io
import os
import copy
INPUT_BASE_FOLDER  = "/home/andrew/Documents/MaterialsDB/pwscf_files/"
RUN_BASE_FOLDER = "/home/andrew/Documents/QE_Runs"
os.makedirs(RUN_BASE_FOLDER,exist_ok=True)

si2_input = "si2.in"

filepath = os.sep.join((INPUT_BASE_FOLDER,si2_input))

run_sub_dir = "Si2_k_sweep"
run_directory = os.sep.join((RUN_BASE_FOLDER,run_sub_dir))
os.makedirs(run_directory,exist_ok=True)

job_init = andrewsalij.io.QEJob(filepath)
job = copy.deepcopy(job_init)
input = job.input

control_dict= input.sections["control"]

structure = input.structure


save_test_path = os.sep.join((run_directory,"si2_saved.in"))
job.save(save_test_path)
