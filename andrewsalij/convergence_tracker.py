

class ConvergenceTester:
    '''Base class for convergence testing some system'''
    def __init__(self,path):
        '''
        Initializes general convergence tester
        :param path: str
        Path to file containing system for convergence testing
        '''
        self.path = path
        self.convergence_list = []
        self.converged_value = None
    def find_convergence(self,convergence_delta):
        '''
        Convergence parameter tolerance in units of parameter
        :param convergence_delta:
        :return:
        '''
        return None

class KPointConvergenceTester(ConvergenceTester):
    '''Subclass for k-point convergence testing'''
    def __init__(self,path):
        super().__init__(path)

