import copy

import numpy as np
import warnings
import os
import andrewsalij.io
'''Convergence tracking of k points'''
class ConvergenceTester:
    '''Base class for convergence testing some system'''
    def __init__(self,path,output_dir = None,job_type = "pwscf"):
        '''
        Initializes general convergence tester
        :param path: str
        :param output_dir: str or None
        :param job_type : str
        Path to file containing system for convergence testing
        '''
        self.path = path
        self.convergence_energy_list = []
        self.convergence_parameter_list = []
        self.output_list = []
        self.converged_energy = None
        self.converged_parameter = None
        if (path):
            self.init_directory, self.filename = os.path.split(self.path)
            self.base_output_directory = self._initialize_convergence_output_dir(output_dir=output_dir)
            self.base_job= andrewsalij.io.load_job(self.path,job_type=job_type)
    def find_convergence(self,convergence_delta):
        '''
        Convergence parameter tolerance in units of parameter
        :param convergence_delta: float
        :return:
        '''
        #TODO write general convergence handling
        return None
    def _initialize_convergence_output_dir(self,output_dir = None):
        '''Initializes the convergence testing by creating output directory. If output_dir = None, makes directory 'conv' in self.path'''
        if output_dir is None:
            output_dir = os.sep.join((self.init_directory,"conv"))
        os.makedirs(output_dir,exist_ok=True)
        return output_dir

class KPointConvergenceTester(ConvergenceTester):
    '''Subclass for k-point convergence testing'''
    def __init__(self,path,output_dir = None,job_type = "pwscf"):
        super().__init__(path,output_dir=output_dir,job_type=job_type)

    @staticmethod
    def create_k_point_array(k_iterator,k_index = 0,effective_weight_array = np.array([1,1,1]),fixed_k_points = None):
        '''
        :param k_iterator:
        :param k_index:
        :param effective_weight:
        :return:
        '''
        k_index,k_iterator = int(k_index), float(k_iterator)
        normalized_weight_array = np.array((effective_weight_array.astype(float))/float(effective_weight_array[k_index])) #normalized to iterator
        k_point_array = (normalized_weight_array*k_iterator).astype(int)
        if fixed_k_points is not None: #for fixed k points, overrides calculated k point array
            try:
                for i in np.arange(3):
                    k_point_force = fixed_k_points[i]
                    if (k_point_force is not None): k_point_array[i] = int(k_point_force)
            except:
                warnings.warn("Manually overriding k points failed: defaulting to initial values")
                pass
        return k_point_array
    def get_effective_k_point_weight(self,weight_type = "uniform"):
        '''
        :param weight_type: str
        :return:
        '''
        if weight_type=="uniform":
            effective_weight_array = np.array([1,1,1]) #no special weighting
        if weight_type=="lattice":
            lattice_array = self.base_job.get_lattice_vectors() #a,b,c (in real space)
            effective_weight_array = 1/lattice_array  #now in reciprocal space
        return effective_weight_array
    def find_convergence(self,convergence_delta,k_iterator_init = int(1),k_step = int(1),k_index = 0,
                         fixed_k_points = None,weight_type = "uniform",run_prefix_str = "",max_iterations = 20,talk = True):
        '''
        Finds the converged energy
        '''
        k_iterator_init,k_step = int(k_iterator_init),int(k_step)
        if (convergence_delta > 0):
            warnings.warn("convergence_delta must be negative. Negativity imposed")
            convergence_delta = -1*np.abs(convergence_delta)
        cur_convergence_delta = convergence_delta-1 # initial value that must be lower that converged value
        cur_k_iterator = k_iterator_init #index for k array construction
        effective_weight_array = self.get_effective_k_point_weight(weight_type = weight_type)
        run_counter = 0 #index for convergence lists
        convergence_flag = False
        while (cur_convergence_delta < convergence_delta): #absolutes taken in case convergence
            if run_counter >= max_iterations:
                print("Max iterations "+str(max_iterations)+" reached. Terminating convergence run")
                break 
            self._run_convergence_iteration(cur_k_iterator,
                                            prefix_str=run_prefix_str,k_index = k_index,
                                            effective_weight_array=effective_weight_array,fixed_k_points=fixed_k_points,talk = talk)
            if (run_counter>=int(1)):
                cur_convergence_delta = self.convergence_energy_list[run_counter]-self.convergence_energy_list[int(run_counter-1)]

            if cur_convergence_delta>0:
                print("Energy increased in final convergence step. Terminating convergence run")
                break
            if (cur_convergence_delta>convergence_delta):
                convergence_flag=True
                break
            cur_k_iterator = cur_k_iterator + k_step
            run_counter = int(run_counter + 1)
        if (convergence_flag):
            self.converged_energy = self.convergence_energy_list[run_counter]
            self.converged_parameter = self.convergence_parameter_list[run_counter]
        else:
            self.converged_energy,self.converged_parameter = None, None
        return self.converged_parameter #returns None if run does not converge
    def _run_convergence_iteration(self,k_iterator,prefix_str = "",k_index = 0,effective_weight_array=np.array([1,1,1]),fixed_k_points = None,talk = True):
        '''Creates directory for convergence iteration, runs jobs, parses output
        :param k_iterator : int
        '''
        iter_k_array = self.create_k_point_array(k_iterator,k_index = k_index,effective_weight_array=effective_weight_array,fixed_k_points=fixed_k_points)
        k_array_str_init  = str(iter_k_array)
        k_array_str = self._clean_k_array_string(k_array_str_init)
        if talk:
            print("Testing K array: "+k_array_str_init)
        iteration_subdir = os.sep.join((self.base_output_directory,"K_"+k_array_str))
        os.makedirs(iteration_subdir,exist_ok=True)
        conv_iter_job = copy.deepcopy(self.base_job)
        conv_iter_job.update_k_points(iter_k_array)
        iteration_output_dir = os.sep.join((iteration_subdir,"out"))
        os.makedirs(iteration_output_dir,exist_ok=True)
        conv_iter_job.update_input_sections({'control':{'outdir':iteration_output_dir}})
        new_input_path = os.sep.join((iteration_subdir,self.filename))
        conv_iter_job.save(new_input_path)
        conv_iter_job.path = new_input_path
        output_filename = conv_iter_job.get_output_filename()
        iteration_output_filepath= os.sep.join((iteration_subdir,output_filename))
        conv_iter_job.run(output_path=iteration_output_filepath,prefix_str=prefix_str)
        iter_output = conv_iter_job.output
        self.output_list.append(iter_output)
        self.convergence_parameter_list.append(iter_k_array)
        self.convergence_energy_list.append(iter_output.final_energy)
        if talk:
            print("Iteration Energy: "+str(iter_output.final_energy))
    @staticmethod
    def _clean_k_array_string(k_array_string,
                              char_list_to_underscore = [","],
                              char_list_to_remove = [" ","(",")","[","]"]):
        '''
        Prepares k_array_string for use in directory names
        :param k_array_string: str
        '''
        k_array_string = k_array_string.lower()
        for char in char_list_to_underscore:
            k_array_string = k_array_string.replace(char,"_")
        for char in char_list_to_remove:
            k_array_string = k_array_string.replace(char,"")
        return k_array_string


