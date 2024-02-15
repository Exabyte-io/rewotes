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
        self.convergence_list = []
        self.converged_value = None
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

        return None
    def _initialize_convergence_output_dir(self,output_dir = None):
        '''Initializes the convergence testing by creating output directory. If output_dir = None, makes directory 'conv' in self.path'''
        if output_dir is None:
            output_dir = os.sep.join((self.init_directory,"conv"))
        os.makedirs(output_dir,exist_ok=True)
        return output_dir

class KPointConvergenceTester(ConvergenceTester):
    '''Subclass for k-point convergence testing'''
    def __init__(self,path):
        super().__init__(path)

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
            lattice_array = self.get_lattice_vectors() #a,b,c (in real space)
            effective_weight_array = 1/lattice_array  #now in reciprocal space
        return effective_weight_array
    def get_lattice_vectors(self):
        '''
        :return:
        '''
