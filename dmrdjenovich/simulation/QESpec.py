
class QESpec:
    """
    Class giving complete info to specify a QuantumEspresso +
    Exabyte + post-analysis run, alongside some helper methods.
    """
    tokens = ["&CONTROL", "&SYSTEM", "&ELECTRONS", "&IONS", "&CELL",
              "&FCP", "&RISM", "ATOMIC_SPECIES", "ATOMIC_POSITIONS",
              "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS",
              "CONSTRAINTS", "ATOMIC_VELOCITIES", "ATOMIC_FORCES",
              "ADDITIONAL_K_POINTS", "SOLVENTS", "HUBBARD"]
    
    def __init__(self, dict):
        """
        Given a dictionary mapping token -> dictionary / list
        This mapped dictionary maps a series of key-value pairs
        The list will be a series of properties
        """
        self.dict = dict
    
    def clone(self):
        """
        Creates a deep copy of the provided QESpec object.
        """
        copy = {}
        for token in QESpec.tokens:
            copy[token] = self.dict[token].copy()
        return QESpec(copy)
        
    def get_lattice(self):
        pass
        
    def set_lattice(self, lattice):
        pass
        
    def get_basis(self):
        pass
        
    def set_basis(self, basis):
        pass
        
    def get_encut(self):
        """
        Returns the planewave energy cutoff of the provided
        QESpec object.
        """
        chunk = self.dict["&SYSTEM"]
        int(chunk.get("ecutwfc"))
        
    def set_encut(self, encut):
        """
        MUTATES this QESpec object with the provided planewave
        energy cutoff.
        """
        chunk = self.dict["&SYSTEM"]
        chunk["ecutwfc"] = encut
        
    def get_k_points(self):
        """
        Assuming automatic k-point generation,
        Returns a list of integers giving the number of k-points
        """
        pass
        
    def set_k_points(self, kpts):
        """
        Assuming automatic k-point generation,
        MUTATES this QESpec with the provided k-point numbers
        """
        pass
        
    @staticmethod
    def parse_file(file):
        """
        Given a filepath (str) containing the simulation input
        parameters, constructs a corresponding QESpec object.
        """
        active = None
        qe_dict = {}
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines:
            
                if len(line) == 0:
                    continue
                if line == "/":
                    continue
                    
                token_line = False
                for t in QESpec.tokens:
                    if t in line:
                        if "&" in t:
                            active = {}
                            qe_dict[t] = active
                        else:
                            active = []
                            qe_dict[t] = active
                            args = line.split(" ")
                            if len(args) > 1:
                                active.append(args[1])
                        token_line = True
                if token_line
                    continue
                
                if "=" in line:
                    parts = line.split("=")
                    active[parts[0]] = parts[1]
                else:
                    active.append(line)
        
        return QESpec(qe_dict)
        
    def write_to_file(self, dest):
        """
        Given a destination filepath (str), writes the appropriate
        input file corresponding to this QESpec.
        """
        with open(dest, "w") as f:
            for token in QESpec.tokens:
                if not token in self.dict:
                    continue
                f.write(token + "\n")
                vals = self.dict[token]
                if "&" in token:
                    f.writelines(map(lambda val : val[0] + "=" + val[1] + "\n", vals.items))
                    f.write("/\n")
                else:
                    f.writelines(vals)
                    f.write("\n")
        
    
                
            
