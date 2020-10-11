
class Material:
    def __init__(self,formula,spacegroup=None,formation_energy=None,E_above_hull=None,
                    volume=None,Nsites=None,density=None,crystal_system=None):
        self.formula = formula
        self.spacegroup = spacegroup
        self.formation_energy = formation_energy
        self.E_above_hull = E_above_hull
        self.volume = volume
        self.Nsites = Nsites
        self.density = density
        self.crystal_system = crystal_system
        self.params = { 
            "spacegroup": self.spacegroup,
            "formation_energy__eV": self.formation_energy, 
            "E_above_hull__eV": self.E_above_hull,
            "volume": self.volume, 
            "Nsites": self.Nsites,
            "density__gm_per_cc": self.density,
            "crystal_system": self.crystal_system, 
        } 
        self.training_params = [ param for (param,value) in self.params.items() if value is not None ]
