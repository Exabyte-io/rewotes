
class Resources:
    """
    Class representing a measure of computational
    work: has time and node / cpu requirements.
    """
    
    def __init__(self, nodes, time):
        self.nodes = nodes
        self.time = time
        
    def get_nodes(self):
        return self.nodes
    
    def get_time(self):
        return self.time
