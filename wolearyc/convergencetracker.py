"""Implementation of convergence tester."""

#Should this really be a package? A python script may be simpler.


#  User-facing interface.
#  User should be able to run the script directly (python ___.py) giving two 
#  arguments:
#     1. Path to input data (a pw.in or POSCAR, INCAR, KPOINTS, POTCAR)
#     2. Optionally, a kinetic energy cutoff (this will override the value in
#         INCAR)
#  
#  The program should print status messages, detailing status of submitted jobs, 
#  verifying job success, giving updates on the convergence, and finally
#  returning the result.
#  
#  """Running jobs."""
#  If I use VASP, the easiest, most flexible option for me (albeit requring lots of
#  user interaction) is to directly copy the VASP input files over to mat3ra, then running
#  the job. Not sure if this is possible (research needed). 
#  Pros of this approach:
#      - complete, explicit control of DFT parameters by the user
#      - great for someone familiar with VASP but not familiar with mat3ra
#  Cons of this approach:
#      - bad for users just starting out with VASP who are not yet familiar with the
#        huge selection of parameters available


"""Class structure."""

def run_vasp_job(input_file_dir, kpoints):
    """ Utility function. Runs a job with vasp.
        Returns...some sort of mat3ra object with all the info?
    """

"""Superclass defining KConverger."""
class KConverger:
    def __init__(self,input_file_dir, initial_kpoints, threshold):
        """ Constructor.

        """
        self.input_file_dir = input_file_dir
        self.initial_kpoints = initial_kpoints
        self.threshold = threshold

    def execute(self):
        """ Executes the convergence test."""

    def is_converged(self,calculation,ref_calculation)
        """ Returns True if less_accurate_calculation is converged based on results of
            more_accurate_calculation. Must be implemented in subclasses.
        """
        raise NotImplementedError

"""Subclass defining KEnergyConverger."""
class KEnergyConverger:
    def __init__(self,input_file_dir,initial_kpoints,threshold):
        super().__init__(input_file_dir,initial_kpoints,threshold)
    
    def is_converged(self,calculation,ref_calculation):
        raise NotImplementedError 


                        