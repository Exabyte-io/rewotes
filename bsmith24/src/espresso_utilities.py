
class Espresso_Calculation:
    """
    This class contains utilities / methods for making and performing calculations with the Quantum Espresso program.
    """

    def __init__(self, template_files_path, ecutwfcs, kpoints):
        self.template_files_path = template_files_path
        self.ecutwfcs = ecutwfcs
        self.kpoints = kpoints

