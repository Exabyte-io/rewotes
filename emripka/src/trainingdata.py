import numpy as np
import json   

def create_id():
    # to-do: need to make this a unique value wrt all user data made and stored
    return "user-%05d" % np.random.randint(0,99999)  

class TrainingData:
    def __init__(self,formula,spacegroup="",formation_energy=0.0,E_above_hull=0.0,band_gap=0.0,has_bandstructure=False,
                    volume=0.0,Nsites=0,theoretical=False,count=0.0,density=0.0,crystal_system=""):
        self.ID = create_id() 
        self.formula = formula
        self.spacegroup = spacegroup 
        self.formation_energy = formation_energy 
        self.E_above_hull = E_above_hull 
        self.band_gap = band_gap 
        self.has_bandstructure = has_bandstructure 
        self.volume = volume 
        self.Nsites = Nsites 
        self.theoretical = theoretical 
        self.count = count 
        self.density = density 
        self.crystal_system = crystal_system 
        self.data_dict = dict() 

    def make_data_dict(self):
        data_dict = {
            self.ID: {
                "formula": self.formula,
                "spacegroup": self.spacegroup,
                "formation_energy__eV": self.formation_energy,
                "E_above_hull__eV": self.E_above_hull,
                "band_gap__eV": self.band_gap,
                "has_bandstructure": self.has_bandstructure,
                "volume": self.volume,
                "Nsites": self.Nsites,
                "theoretical": self.theoretical,
                "count": self.count,
                "density__gm_per_cc": self.density,
                "crystal_system": self.crystal_system 
            }
        }
        self.data_dict = data_dict

    def show_data(self):
        self.make_data_dict()
        print(self.data_dict)

    def store_data(self):
        """ This method will store the user training data to the data directory."""
        self.make_data_dict()

        json_path = "../data/training/materialsproject_json/"

        # open the existing user_data.json file
        with open(f"{json_path}/user_data.json","r") as fname:
            try:
                user_data = dict(json.load(fname))
            except:
                user_data = dict() 
        fname.close()

        user_data[self.ID] = self.data_dict[self.ID]
        with open(f"{json_path}/user_data.json","w") as fname:
            json.dump(user_data,fname,indent=4)
        fname.close()
