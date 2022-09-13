import xml.etree.ElementTree as ET

class QEAnalysis(object):
    """
    Class containing standard methods for processing the output
    of a Quantum Espresso simulation.
    """
    
    def __init__(self, src):
        self.root = ET.parse(src).getroot()
        
    def get_energy(self):
        """
        Returns an array of energies, each corresponding to an
        ionic relaxation step in order of initial to final
        geometry.
        """
        raw = self.root.findall("output/total_energy/etot")
        return [float(x.text) for x in raw]

