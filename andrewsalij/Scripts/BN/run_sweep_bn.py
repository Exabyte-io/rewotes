import os
import andrewsalij.convergence_tracker as convergence_tracker

'''Script for testing that k sweeping is working for boron nitride.'''

INPUT_BASE_FOLDER  = "/home/andrew/Documents/MaterialsDB/pwscf_files/"
RUN_BASE_FOLDER = "/home/andrew/Documents/QE_Runs"
os.makedirs(RUN_BASE_FOLDER,exist_ok=True)

compound_str = "br4ca2"

input = compound_str+".in"
filepath = os.sep.join((INPUT_BASE_FOLDER,input))

run_sub_dir = compound_str+"_k_sweep"
run_directory = os.sep.join((RUN_BASE_FOLDER,run_sub_dir))
convergence_tester = convergence_tracker.KPointConvergenceTester(filepath,output_dir = run_directory)

k_array_conv = convergence_tester.find_convergence(-.001,run_prefix_str="mpirun -np 8",
                                                   k_index=2,k_step=2,k_iterator_init=3)

convergence_tester.make_report_figure(compound_str+"_k_converged.png",x_axis_type="convergence_parameter")