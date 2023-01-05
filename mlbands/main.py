from mp_api.client import MPRester
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# from mendeleev import element
# elements = {}
# for i in range(1,119):
# 	elements[element(i).symbol] = i
elements={'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}
 
from mlbands.misc import *

import numpy as np
import matplotlib.pyplot as plt

class Material:
    def __init__( self, api_key, material_ID = 'mp-1103503', box_array = None):
        '''Material class. To store and process material data

        Parameters
        ----------
        API_KEY : str
            secret API KEY used to run MPRester data requests
        material_ID : int
            ID of selected material (without "mp-" prefix)
        box_array : (float, float, float)
            3-D Tensor representing material chemical env. 
        '''
        self.API_KEY = api_key
        self.material_ID = material_ID
        self.box_array = box_array          # material already transformed to box data form (optional)

    def bands(self,nonzero_gap=False):
        '''obtain band gap information for requested Material 

        Parameters
        ----------
        nonzero_gap: bool
            if True, only returns nonzero band gaps 
        '''
        with MPRester(api_key=self.API_KEY) as mpr:
            #adapted from https://matsci.org/t/obtain-large-numbers-of-band-structures/3780
            bandstructure = None 
            try:
                bandstructure = mpr.get_bandstructure_by_material_id(self.material_ID,line_mode=False)
            except:
                pass
            if bandstructure:
                band_gap = bandstructure.get_band_gap()
                print('Band Gap: {} eV\nDirect Gap: {}\nMetallic: {}'.\
                    format(band_gap['energy'],\
                        'Yes' if band_gap['direct'] else 'No',\
                        'No' if band_gap['transition'] else 'Yes'))
                        
                if nonzero_gap:
                    if band_gap['energy']: return band_gap
                    else: return 0
                else: return band_gap
            else:
                return 0

    def load_structure(self, conventional=True):
        '''helper function to load structural information of Material

        Parameters
        ----------
        conventional: bool
            if True, returns conventional standard structure else returns primitive
        '''
        with MPRester(self.API_KEY) as mpr:

            # first retrieve the relevant structure
            structure = mpr.get_structure_by_material_id(self.material_ID)
        
        # important to use the conventional structure to ensure
        # that peaks are labelled with the conventional Miller indices
        sga = SpacegroupAnalyzer(structure)
        conventional_structure = sga.get_conventional_standard_structure()
        
        if conventional:
            structure = conventional_structure

        return structure

       
    def structural(self):
        '''prints structural information
        '''
        structure = Material(self.API_KEY, self.material_ID).load_structure()

        print(structure.lattice)
        print(structure.sites)  #https://pymatgen.org/pymatgen.core.sites.html?highlight=periodicsite#pymatgen.core.sites.PeriodicSite

        Nsites = len(structure.sites)
        for i in range(Nsites):
            print('\n\n')
            # print(structure.sites[i])
            print(structure.sites[i].species)
            print(structure.sites[i].coords)
            print(structure.sites[i].frac_coords)     
    

    def to_xyz(self, fractional=False):
        '''processes the material structure as an array with an xyz coordinate format 
        i.e. [atom_numbers, x_coord, y_coord, z_coord]

        Parameters
        ----------
        fractional: bool
            if False, returns xyz coordinates else returns fractional (a,b,c) coordinates
        '''
        structure = Material(self.API_KEY, self.material_ID).load_structure(conventional=True)

        Nsites = len(structure.sites)

        xyz_array = np.zeros((Nsites,4))

        if not fractional:
            for i in range(Nsites):
                # atom_num = elements(str(structure.sites[i].specie)).atomic_number # specie is not a typo
                atom_num = elements[str(structure.sites[i].specie)] # specie is not a typo

                xyz_array[i] = [atom_num,*structure.sites[i].coords]
        else:
            for i in range(Nsites):
                # atom_num = elements(str(structure.sites[i].specie)).atomic_number 
                atom_num = elements[str(structure.sites[i].specie)] 

                xyz_array[i] = [atom_num,*structure.sites[i].frac_coords]

        return xyz_array


    def to_box(self, fractional=False):
        '''processes the material structure as a sparse, 3-D Tensor "Box" with atom_numbers as the 
        non-zero entries
        i.e. box = Tensor((z_coord, y_coord, x_coord))

        Parameters
        ----------
        fractional: bool
            if False, returns tensor space from xyz coordinates else from fractional (a,b,c) coordinates
        '''
        xyz_array = Material(self.API_KEY, self.material_ID).to_xyz(fractional)

        coords = xyz_array[:,1:]            
        MAX = np.ceil(np.max(coords)).astype('int')
        MIN = np.floor(np.min(coords)).astype('int')
        
        xyz_array[:,1:] -= MIN

        L = (MAX-MIN+1)

        box = np.zeros((L,L,L))
        for i in xyz_array:
            print(i)

            atom,x,y,z = i.astype('int')
            box[z,y,x] = atom

        return box
            
    def visual(self, spacing = 1, fractional=False):
        '''simple visualization tool for Material (plot)

        Parameters
        ----------
        spacing: float
            increases spacing by multiplying coordinate distance
        fractional: bool
            if False, returns tensor space from xyz coordinates else from fractional (a,b,c) coordinates
        '''
        ax = plt.axes(projection='3d')
        colors = np.linspace(2**20,2**24,118,dtype='int') #divide color range into 118 colors (for the 118 chemical elements)


        if self.box_array is not None:
            #presence of box_array supplants material_ID
            for i in np.argwhere(self.box_array):
                x,y,z = i
                atom = int(self.box_array[tuple(i)])
                ax.scatter3D(x,y,z, s=100, c="#"+hex(colors[atom])[2:])

        else:
            xyz_array = Material(self.API_KEY, self.material_ID).to_xyz(fractional)
            xyz_array[:,1:]*=spacing

            for i in xyz_array:
                atom,*xyz = i.astype('int')
                ax.scatter3D(*xyz, s=100, c="#"+hex(colors[atom])[2:])


        set_axes_equal(ax)           

        plt.axis('off')
        plt.show()


    def XRD(self):
        '''XRD pattern data for Material
        '''
        structure = Material(self.API_KEY, self.material_ID).load_structure(self.API_KEY)
        
        # this example shows how to obtain an XRD diffraction pattern
        # these patterns are calculated on-the-fly from the structure
        calculator = XRDCalculator(wavelength='CuKa')
        pattern = calculator.get_pattern(structure)

        print('\npattern:\n',pattern)
        

    def thermo(self):   
        '''thermodynamic data for Material
        '''
        with MPRester(api_key=self.API_KEY) as mpr:

            # for a single material
            # thermo_docs = mpr.thermo.get_data_by_id(self.material_ID)      #   DOES NOT WORK msg(Item with thermo_id = mp-1103503 not found)
            
            # for many materials, it's much faster to use
            # the `search` method, where additional material_ids can 
            # be added to this list
            thermo_docs = mpr.thermo.search(material_ids=[self.material_ID])

        # print(thermo_docs[0].energy_per_atom)
        return thermo_docs

    def magnetism(self):      
        '''magnetic properties data for Material
        '''
        with MPRester(api_key=self.API_KEY) as mpr:
            magnetism_doc = mpr.magnetism.get_data_by_id('mp-1103503')

        # print(magnetism_doc.total_magnetization)
        return magnetism_doc

class Group:
    def __init__(self,api_key):
        '''Group class for bulk calculations

        Parameters
        ----------
        API_KEY : str
            secret API KEY used to run MPRester data requests
        materials : int
            IDs of selected material (with "mp-" prefix)
        X : [varies] (float,float,float); ((float,float,float), *float); (*float)
            quantitative predictors ("X-value")
        Y : float
            band gap (eV) quantitative response ("Y-value","label")
        '''
        self.API_KEY = api_key
        self.materials = []         # material IDs (materials found in ID_list; see data_make inputs)
        self.X = []                 # quantitative predictors ("X-value")
        self.Y = []                 # quantitative response ("Y-value","label")
        
        self.box_lengths = []       # length scale of boxes (3-D Tensors of material chemical env.)
        self.max_length = 0

    def transfer(self, loaded_data):
        '''transfer constructor variable information to Group class from loaded data

        Parameters
        ----------
        loaded data : <"Group" class object>
            file containing processed Materials in a Group container
        '''
        self.X = loaded_data.X
        self.Y = loaded_data.Y
        self.materials = loaded_data.materials

    def data_make(self, ID_list = range(1,10), nonzero_gap=False ):    
        '''generate Group data for selected ID_list numbers [ for materials which exist with such IDs ]

        Parameters
        ----------
        ID_list : list[int]
            list of material_ID numbers (i.e. without "mp-" prefix)
        nonzero_gap: bool
            if True, only returns nonzero band gaps 
        '''
        for i in ID_list:
            material = Material(self.API_KEY, 'mp-'+str(i))
            BG = material.bands(nonzero_gap)
            if BG:                              # if material with material_ID exists
                self.materials.append('mp-'+str(i))

                # quantitative response ("Y-value","label")
                self.Y.append(BG['energy'])     # append band gap energy ( eV ) 

                box = material.to_box()
                self.box_lengths.append(box.shape[0])
                # quantitative predictors ("X-value")
                self.X.append(box)              # append 3-D Tensor representing material chemical env. 
                # for extra_x in features:
                #     self.X.append(extra_x)

        self.Y = np.array(self.Y)

    def data_expand(self,boxes=True,*property_funcs):
        '''expand "X-values" data to other chemical characteristics 

        e.g. data_expand(True, thermo, magnetic) appends "thermo" and "magnetic" characteristics to X in addition to 3-D Tensors
        e.g. data_expand(False, thermo, magnetic) creates a new X-value object and appends "thermo" and "magnetic" characteristics alone

        Parameters
        ----------
        boxes : bool
            removes 3-D tensors to X variable if False, 
        *property_funcs: <function_1>, ..., <function_N>
            operates functions on Materials which extract additional properties and appends them to X
        '''
        boxes = self.X
        self.X = [] 
        # load a property for all materials
        for i in self.materials:
            material = Material(self.API_KEY, 'mp-'+str(i))

            k=0
            material_props = []
            material_props.append(boxes[k]) if boxes else None
            for func in property_funcs:
                material_props.append(material.func)
            
            self.X.append(material_props)
            k+=1

       
    def resize_boxes(self, L=32):
        '''resize 3-D tensor boxes to fit neural network architecture (32x32x32)

        Parameters
        ----------
        L: int
            length-scale of boxes (3-D Tensors representing material chemical env.)
        '''
        max_length = np.max(self.box_lengths)

        if L >= max_length:


            self.X = np.array([  np.pad(self.X[i],\
                        ( (0,L-self.box_lengths[i]),(0,L-self.box_lengths[i]),(0,L-self.box_lengths[i]) ) 
                        ) for i in range(len(self.X)) ])

        # print(self.X.shape)

    


