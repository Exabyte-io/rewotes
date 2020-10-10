
class Material:
    def __init__(self,formula,formation_energy=None,volume=None,density=None):
        self.formula = formula
        self.formation_energy = formation_energy
        self.volume = volume
        self.density = density
        self.params = { 
            "density__gm_per_cc": self.density,
            "volume": self.volume, 
            "formation_energy_eV": self.formation_energy, 
        } 
        self.training_params = [ param for (param,value) in self.params.items() if value is not None ]
