from pymatgen.io import pwscf

'''Input and output handling
Formats supported: Quantum Espresso (PWscf input)
'''
#TODO: add additional format support

class Job():
    '''
    Base class for a file to run calculations
    '''
    def __init__(self,path,to_load = False):
        self.path = path
class QEJob(Job):
    '''
    Class for running Quantum Espresso calculations. Extends Job()
    '''
    def __init__(self,path,to_load = False):
        super().__init__(path,to_load)
        self.input = self.load(path) #PWInput object
    @staticmethod
    def load(path):
        return pwscf.PWInput.from_file(path)

