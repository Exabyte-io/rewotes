import math
import numpy as np

class QESpec(object):
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
        
    def get_lattice_type(self):
        """
        Returns the type of lattice specification used for this
        simulation, an integer code.
        """
        return int(self.dict["&SYSTEM"].get("ibrav", 0))
        
    def get_lattice(self):
        """
        Returns the generators of the lattice used for this simulation
        as column vectors of 2-dimensional array.
        Units are in angstroms ALWAYS
        """
        BOHR_CONVERSION = 0.529177
        cell_type = int(self.dict["&SYSTEM"].get("ibrav", 0))
        cell_dim_1 = self.dict["&SYSTEM"].get("celldm(1)", None)
        
        if cell_type == 0:
            SRC = self.dict.get("CELL_PARAMETERS", None)
            if SRC is None:
                raise Exception("Unable to parse lattice.")
            if cell_dim_1 == "A":
                factor = 1
            elif cell_dim_1 is not None:
                factor = BOHR_CONVERSION*float(cell_dim_1)
            else:
                factor = BOHR_CONVERSION
            if SRC[0] == "bohr":
                factor = BOHR_CONVERSION
            elif SRC[0] == "angstrom":
                factor = 1
            generators = []
            generators.append([float(x)*factor for x in SRC[1].strip().split(" ")])
            generators.append([float(x)*factor for x in SRC[2].strip().split(" ")])
            generators.append([float(x)*factor for x in SRC[3].strip().split(" ")])
            return np.transpose(generators)
            
        if cell_dim_1 is None or cell_dim_1 == "A":
            A = float(self.dict["&SYSTEM"]["A"])
            B = float(self.dict["&SYSTEM"]["B"])
            C = float(self.dict["&SYSTEM"]["C"])
            cell_dim_1 = A
            cell_dim_2 = B/A
            cell_dim_3 = C/A
            cell_dim_4 = float(self.dict["&SYSTEM"]["cosAB"])
            cell_dim_5 = float(self.dict["&SYSTEM"]["cosAC"])
            cell_dim_6 = float(self.dict["&SYSTEM"]["cosBC"])
        elif cell_type != 0:
            cell_dim_1 = float(cell_dim_1)*BOHR_CONVERSION
            cell_dim_2 = float(self.dict["&SYSTEM"].get("celldm(2)", -1))
            cell_dim_3 = float(self.dict["&SYSTEM"].get("celldm(3)", -1))
            cell_dim_4 = float(self.dict["&SYSTEM"].get("celldm(4)", -1))
            cell_dim_5 = float(self.dict["&SYSTEM"].get("celldm(5)", -1))
            cell_dim_6 = float(self.dict["&SYSTEM"].get("celldm(6)", -1))
        else:
            raise Exception("Unable to parse lattice")
            
        if cell_type == 1:
            # SIMPLE CUBIC
            a = cell_dim_1
            return [[a, 0, 0],[0, a, 0],[0, 0, a]]
        if cell_type == 2:
            # FACE-CENTERED-CUBIC
            a = cell_dim_1 / 2
            return [[-a, 0, -a],[0, a, a],[a, a, 0]]
        if cell_type == 3:
            # BODY-CENTERED-CUBIC
            a = cell_dim_1 / 2
            return [[a, -a, -a],[a, a, -a],[a, a, a]]
        if cell_type == -3:
            # BODY-CENTERED-CUBIC
            a = cell_dim_1 / 2
            return [[-a, a, a],[a, -a, a],[a, a, -a]]
        if cell_type == 4:
            # PRIMITIVE-HEXAGONAL-TRIGONAL
            a = cell_dim_1
            r = cell_dim_3
            return [[a, -a/2, 0],[0, a*sqrt(3)/2, 0],[0, 0, a*r]]
        if cell_type == 5:
            # RHOMBOHEDRAL-TRIGONAL
            a = cell_dim_1
            c = cell_dim_4
            tx = sqrt((1-c)/2)
            ty = sqrt((1-c)/6)
            tz = sqrt((1+2*c)/3)
            return [[a*tx, 0, -a*tx],[-a*ty, 2*a*ty, -a*ty],[a*tz, a*tz, a*tz]]
        if cell_type == -5:
            # RHOMBOHEDRAL-TRIGONAL
            a = cell_dim_1 / sqrt(3)
            c = cell_dim_4
            tx = sqrt((1-c)/2)
            ty = sqrt((1-c)/6)
            tz = sqrt((1+2*c)/3)
            u = tz - 2*sqrt(2)*ty
            v = tz + sqrt(2)*ty
            return [[a*u, a*v, a*v], [a*v, a*u, a*v], [a*v, a*v, a*u]]
        if cell_type == 6:
            # PRIMITIVE-TETRAGONAL
            a = cell_dim_1
            r = cell_dim_3
            return [[a, 0, 0], [0, a, 0], [0, 0, a*r]]
        if cell_type == 7:
            # BODY-CENTERED-TETRAGONAL
            a = cell_dim_1 / 2
            r = cell_dim_3
            return [[a, a, -a], [-a, a, -a], [a*r, a*r, a*r]]
        if cell_type == 8:
            # PRIMITIVE-ORTHORHOMBIC
            a = cell_dim_1
            r1 = cell_dim_2
            r2 = cell_dim_3
            return [[a, 0, 0], [0, a*r1, 0], [0, 0, a*r2]]
        if cell_type == 9:
            # BASE-CENTERED-ORTHORHOMBIC
            a = cell_dim_1
            r1 = cell_dim_2
            r2 = cell_dim_3
            return [[a/2, -a/2, 0], [a*r1/2, a*r1/2, 0], [0, 0, a*r2]]
        if cell_type == -9:
            # BASE-CENTERED-ORTHORHOMBIC
            a = cell_dim_1
            r1 = cell_dim_2
            r2 = cell_dim_3
            return [[a/2, a/2, 0], [-a*r1/2, a*r1/2, 0], [0, 0, a*r2]]
        if cell_type == 91:
            # BASE-CENTERED-ORTHORHOMBIC
            a = cell_dim_1
            r1 = cell_dim_2
            r2 = cell_dim_3
            return [[a, 0, 0], [0, a*r1/2, a*r1/2], [0, -a*r2/2, a*r2/2]]
        if cell_type == 10:
            # FACE-CENTERED-ORTHORHOMBIC
            a = cell_dim_1
            r1 = cell_dim_2
            r2 = cell_dim_3
            return [[a/2, a/2, 0], [0, a*r1/2, a*r1/2], [a*r2/2, 0, a*r2/2]]
        if cell_type == 11:
            # BODY-CENTERED-ORTHORHOMBIC
            a = cell_dim_1
            r1 = cell_dim_2
            r2 = cell_dim_3
            return [[a/2, a*r1/2, a*r2/2], [-a/2, a*r1/2, a*r2/2], [-a/2, -a*r1/2, a*r2/2]]
        if cell_type == 12:
            # PRIMITIVE-MONOCLINIC
            a = cell_dim_1
            r1 = cell_dim_2
            r2 = cell_dim_3
            angle = math.acos(cell_dim_4)
            return [[a, a*r1*math.cos(angle), 0], [0, a*r1*math.sin(angle), 0], [0, 0, a*r2]]
        if cell_type == -12:
            # PRIMITIVE-MONOCLINIC
            a = cell_dim_1
            r1 = cell_dim_2
            r2 = cell_dim_3
            angle = math.acos(cell_dim_5)
            return [[a, 0, a*r2*math.cos(angle)], [0, a*r1, 0], [0, 0, a*r2*math.sin(angle)]]
        if cell_type == 13:
            # BASE-CENTERED-MONOCLINIC
            a = cell_dim_1
            r1 = cell_dim_2
            r2 = cell_dim_3
            angle = math.acos(cell_dim_4)
            return [[a/2, a*r1*math.cos(angle), a/2], [0, a*r1*math.sin(gamma), 0], [-a*r2/2, 0, a*r2/2]]
        if cell_type == -13:
            # BASE-CENTERED-MONOCLINIC
            a = cell_dim_1
            r1 = cell_dim_2
            r2 = cell_dim_3
            angle = math.acos(cell_dim_5)
            return [[a/2, -a/2, a*r2*math.cos(angle)], [a*r1/2, a*r1/2, 0], [0, 0, a*r2*math.sin(angle)]]
        if cell_type == 14:
            # TRICLINIC
            a = cell_dim_1
            r1 = cell_dim_2
            r2 = cell_dim_3
            alpha = math.acos(cell_dim_4)
            beta = math.acos(cell_dim_5)
            gamma = math.acos(cell_dim_6)
            return [[a, a*r1*math.cos(gamma), a*r2*math.cos(beta)],
                    [0, a*r1*math.sin(gamma), a*r2*(math.cos(alpha) - math.cos(beta)*math.cos(gamma))/math.sin(gamma)],
                    [0, 0, a*r2*sqrt(1 + 2*math.cos(alpha)*math.cos(beta)*math.cos(gamma) - math.cos(alpha**2 - math.cos(beta)**2 - math.cos(gamma)**2)/math.sin(gamma))]]
        
    def set_lattice(self, generators):
        """
        Sets the generators of the lattice used for this simulation
        passed in as column vectors of a 2-dimensional array.
        Units are in angstroms ALWAYS
        """
        self.dict["&SYSTEM"]["ibrav"] = "0"
        self.dict["&SYSTEM"]["celldm(1)"] = "A"
        build = []
        self.dict["CELL_PARAMETERS"] = build
        build[0] = "angstrom"
        for i in range(0, len(generators)):
            for j in range(0, len(generators)):
                build[i] += str(generators[j][i]) + " "
        
    def get_species(self):
        pass
        
    def set_species(self, species):
        pass
        
    def get_positions(self):
        """
        Returns the atomic motif in relative coordinates
        of the crystal [[species],[positions]].
        Positions are arranged as column vectors of a
        matrix.
        """
        SRC = self.dict["ATOMIC_POSITIONS"]
        if SRC[0] == "alat":
            raise Exception("Not yet implemented.")
        if SRC[0] == "crystal_sg":
            raise Exception("Not yet implemented.")
        
        species = []
        loc = []
        for i in range(1, len(SRC)):
            tokens = SRC[i].strip().split(" ")
            species.append(tokens[0])
            loc.append([float(x) for x in tokens[1:]])
        loc = list(np.transpose(loc))
        for i in range(0, len(loc)):
            loc[i] = list(loc[i])
        if SRC[0] == "crystal":
            return [species, loc]
        
        mult = np.linalg.inv(self.get_lattice())
        loc = np.matmul(mult, loc)
        loc = list(loc)
        for i in range(0, len(loc)):
            loc[i] = list(loc[i])
        if SRC[0] == "angstrom":
            return [species, loc]
        BOHR_CONVERSION = 0.529177
        if SRC[0] == "bohr":
            for i in range(0, len(loc)):
                for j in range(0, len(loc[i])):
                    loc[i][j] *= BOHR_CONVERSION
            return [species, loc]
        
    def set_positions(self, species, positions):
        """
        Sets the atomic motif in relative coordinates
        of the crystal, positions are specified as
        column vectors.
        """
        SRC = []
        self.dict["ATOMIC_POSITIONS"] = SRC
        SRC[0] = "crystal"
        transp = np.transpose(positions)
        for i in range(0, len(species)):
            SRC[i+1] = species[i] + " " + " ".join((str(transp[i][x]) for x in range(0, len(transp[i]))))
        
    def get_encut(self):
        """
        Returns the planewave energy cutoff of the provided
        QESpec object.
        """
        chunk = self.dict["&SYSTEM"]
        float(chunk.get("ecutwfc"))
        
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
        SRC = self.dict["K_POINTS"]
        if SRC[0] != "automatic":
            raise Exception("Parsing of these kpoints not yet implemented")
        return [int(x) for x in SRC[1].strip().split(" ")]
        
    def set_k_points(self, kpts):
        """
        Assuming automatic k-point generation,
        MUTATES this QESpec with the provided k-point numbers
        """
        build = []
        self.dict["K_POINTS"] = build
        build[0] = "automatic"
        build[1] = ""
        for kpt in kpts:
            build[1] += str(kpt) + " "
        build[1].trim()
        
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
                if token_line:
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
                if not token in self.dict.keys():
                    continue
                vals = self.dict[token]
                if "&" in token:
                    f.write(token + "\n")
                    f.writelines(map(lambda val : val[0] + "=" + val[1] + "\n", vals.items))
                    f.write("/\n")
                else:
                    if " " not in vals[0]:
                        f.write(token + " " + vals[0] + "\n")
                        f.writelines(vals[1:])
                        f.write("\n")
                    else:
                        f.write(token + "\n")
                        f.writelines(vals)
                        f.write("\n")
        
    
                
            
