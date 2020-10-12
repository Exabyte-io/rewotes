class Calculation(object):
    def __init__(self, incar: str, poscar: str, kpoints: str = None):
        """
        Object containing methods related to a VASP calculation.

        :param incar: Path to the INCAR
        :type incar: str
        :param poscar: Path to the POSCAR
        :type poscar: str
        :param kpoints: Path to the KPOINTS
        :type kpoints: str
        """
        self.incar = incar
        self.poscar = poscar
        self.kpoints = kpoints

def example_function():
    """
    This is a test
    :return:
    """
    return None

if __name__ == "__main__":
    pass