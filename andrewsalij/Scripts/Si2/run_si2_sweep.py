import andrewsalij.io
import os
import copy
import andrewsalij.convergence_tracker as convergence_tracker


'''Script for testing that Si2 k sweeping is working'''

INPUT_BASE_FOLDER  = "/home/andrew/Documents/MaterialsDB/pwscf_files/"
RUN_BASE_FOLDER = "/home/andrew/Documents/QE_Runs"
os.makedirs(RUN_BASE_FOLDER,exist_ok=True)

si2_input = "si2.in"


filepath = os.sep.join((INPUT_BASE_FOLDER,si2_input))

run_sub_dir = "Si2_k_sweep"
run_directory = os.sep.join((RUN_BASE_FOLDER,run_sub_dir))
convergence_tester = convergence_tracker.KPointConvergenceTester(filepath,output_dir = run_directory)

base_job = convergence_tester.base_job
update_job = copy.deepcopy(base_job)
sections = update_job.input.sections


convergence_tester.find_convergence(-.005)

convergence_list = convergence_tester.convergence_parameter_list
convergence_energy_list = convergence_tester.convergence_energy_list

convergence_tester.make_report_figure("si2_k_converged.png",x_axis_type="convergence_parameter")

