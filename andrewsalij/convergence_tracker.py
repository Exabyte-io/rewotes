import copy
import matplotlib.pyplot as plt
import numpy as np
import warnings
import os
import andrewsalij.io
'''Convergence tracking of parameters, particularly testing for k point convergence.'''
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
        '''Initializes KPointConvergenceTester as an instance of ConvergenceTester
        :param path: str
        :param output_dir: str or None
        :param job_type: str
            "pwscf": Quantum Espresso calculation (pw.x)
        '''
        super().__init__(path,output_dir=output_dir,job_type=job_type)

    @staticmethod
    def create_k_point_array(k_iterator,k_index = 0,effective_weight_array = np.array([1,1,1]),fixed_k_points = None):
        '''
        Static method for the creation of a np.ndarray (size 3). See find_convergence() for details
        :param k_iterator: int
        :param k_index: int
        :param effective_weight_array: np.ndarray (size 3) (default: np.array([1,1,1]))
            Array for weighting of k indices. Defaults to uniform weighting.
        :param fixed_k_points: (value, value, value) where value is int or None
        :return: np.ndarray (size 3, type int)
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
    def _get_effective_k_point_weight(self,weight_type = "uniform"):
        '''
        See find_convergence() for weight_type descriptions
        :param weight_type: str
        :return: np.ndarray (size 3, type int)
        '''
        if weight_type=="uniform":
            effective_weight_array = np.array([1,1,1]) #no special weighting
        if weight_type=="lattice":
            lattice_array = self.base_job.get_lattice_vectors() #a,b,c (in real space)
            effective_weight_array = 1/lattice_array  #now in reciprocal space
        return effective_weight_array
    def find_convergence(self,convergence_delta,k_iterator_init = int(1),k_step = int(1),k_index = 0,
                         fixed_k_points = None,weight_type = "uniform",run_prefix_str = "",max_iterations = 20,talk = True,
                         force_energy_decrease = False):
        '''
        Finds the converged energy for a given convergence delta for the input from the initialization.
        Creates and runs in subdirectories within self.output_dir scf calculations for target system.
        Convergence defined where the latest decrease in energy is less than the provided convergence_delta.
        Convergence iterates over a single k point index (default: first lattice dimension) from k_iterator_init (default: 1)
        in increments of k_step (default: 1). Terminates after set max_iterations (default: 20).
        :param convergence_delta : float (in units of method energy)
            Tolerance for convergence
        :param k_iterator_init : int (default: 1)
            Initial value to begin the iterated k_index
        :param k_step: int  (default: 1)
            Value to increase the k_iterator by each iteration
        :param k_index: int (default: 0)
            Scalar index of k points to iterate over
        :param fixed_k_points: (value, value, value) where value is int or None
            Optional parameter to fix certain k point values to arbitrary integers.
            None either as a whole parameter or in the given tuple index implies that k point
            value increases with iteration.
        :param weight_type: str (default: "uniform")
            Scaling type for k point iteration:
            "uniform":
                All free (not in fixed_k_points) indices are identical
            "lattice":
                K points scaled according to integers closed to reciprocal lattice vectors
        :param run_prefix_str (default: "")
            String to inset (e.g. "mpirun -np 4", before method call (e.g., "pw.x")).
            Enables support for method parallelism
        :param max_iterations: int (default: 20)
            Number of iterations to run before killing convergence
        :param talk: bool (default: True)
            If True, outputs progress to console
        :param force_energy_decrease: bool (default: False)
            If True, kills convergence if energy ever increases
        :return: np.ndarray (size 3, type int) or None
            returns None if run failed to converge
        '''
        k_iterator_init,k_step = int(k_iterator_init),int(k_step)
        if (convergence_delta > 0):
            warnings.warn("convergence_delta must be negative. Negativity imposed")
            convergence_delta = -1*np.abs(convergence_delta)
        cur_convergence_delta = convergence_delta-1 # initial value that must be lower that converged value
        cur_k_iterator = k_iterator_init #index for k array construction
        effective_weight_array = self._get_effective_k_point_weight(weight_type = weight_type)
        run_counter = 0 #index for convergence lists
        convergence_flag = False
        while (not convergence_flag): #absolutes taken in case convergence
            if run_counter >= max_iterations:
                print("Max iterations "+str(max_iterations)+" reached. Terminating convergence run")
                break 
            self._run_convergence_iteration(cur_k_iterator,
                                            prefix_str=run_prefix_str,k_index = k_index,
                                            effective_weight_array=effective_weight_array,fixed_k_points=fixed_k_points,talk = talk)
            if (run_counter>=int(1)):
                cur_convergence_delta = self.convergence_energy_list[run_counter]-self.convergence_energy_list[int(run_counter-1)]

            if cur_convergence_delta>0 and force_energy_decrease:
                print("Energy increased in final convergence step. Terminating convergence run")
                break
            if (cur_convergence_delta>convergence_delta and cur_convergence_delta<0):
                print("Run converged")
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
    def make_report_figure(self,filename,output_directory = None,x_axis_type = "iteration_number",x_label = "Convergence Iteration",y_label = "Energy (Ry)",
                           plot_params = {},to_show = True):
        '''
        Creates basic report figure for k point convergence testing. Plot parameters may be
        overridden via matplotlib.pyplot.rcParams.update()
        '''
        fig, ax = plt.subplots()
        if output_directory is not None:
            full_filename = os.sep.join((output_directory,filename))
        else:full_filename = filename
        n_its = len(self.convergence_energy_list)
        iteration_array = np.linspace(1, n_its, n_its)
        energy_array = np.array(self.convergence_energy_list)
        ax.plot(iteration_array,energy_array,**plot_params)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        if (x_axis_type == "iteration_number"):
            #No additional functionality at present--kept for case switching
            dummy_var = None
        elif (x_axis_type == "convergence_parameter"):
            x_ticks = ax.get_xticks()[1:-1]# removing edge indices that give left and right bounds of axis
            x_labels = self.convergence_parameter_list
            x_tick_indices = np.unique((x_ticks-1).astype(int)) #shifting index from iteration 1 to pythonic 0
            x_ticks = x_tick_indices+1 #resetting after removing duplicates
            x_labels_subset = [x_labels[idx] for idx in list(x_tick_indices)]
            ax.set_xticks(x_ticks,labels = x_labels_subset)
        else:
            ValueError("Invalid x_axis_type: "+str(x_axis_type))
        plt.tight_layout()
        fig.savefig(full_filename)
        if to_show: fig.show()


