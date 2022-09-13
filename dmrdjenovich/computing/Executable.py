
class Executable(object):
    """
    Class representing a chunk of computational work.
    
    Has time and node / cpu requirements, specified
    by a Resources object.
    """
    
    def get_resources(self):
        """
        Returns the required resources for this
        computational task.
        """
        pass
        
    def run(self, envr):
        """
        Given a Resource object representing the
        available computational resources, runs
        this computational task.
        """
        pass
