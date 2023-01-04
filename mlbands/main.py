from mp_api.client import MPRester
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from mendeleev import element

from mlbands.misc import *

import numpy as np
import matplotlib.pyplot as plt

class Material:
    def __init__(self):
        self.API_KEY = ''
        self.structure_ID = 'mp-1103503'

    def bands(self):
        with MPRester(api_key=self.API_KEY) as mpr:
            #adapted from https://matsci.org/t/obtain-large-numbers-of-band-structures/3780
            bandstructure = None 
            try:
                bandstructure = mpr.get_bandstructure_by_material_id(self.structure_ID,line_mode=False)
            except:
                pass
            if bandstructure:
                band_gap = bandstructure.get_band_gap()
                
                print('Band Gap: {} eV\nDirect Gap: {}\nMetallic: {}'.\
                    format(band_gap['energy'],\
                        'Yes' if band_gap['direct'] else 'No',\
                        'No' if band_gap['transition'] else 'Yes'))


    def load_structure(self, api_key, conventional=True):
        with MPRester(api_key) as mpr:
            # first retrieve the relevant structure
            structure = mpr.get_structure_by_material_id(self.structure_ID)
        
        # important to use the conventional structure to ensure
        # that peaks are labelled with the conventional Miller indices
        sga = SpacegroupAnalyzer(structure)
        conventional_structure = sga.get_conventional_standard_structure()
        
        if conventional:
            structure = conventional_structure

        return structure

       
    def structural(self):

        structure = Material().load_structure(self.API_KEY)


        print(structure.lattice)
        print(structure.sites)  #https://pymatgen.org/pymatgen.core.sites.html?highlight=periodicsite#pymatgen.core.sites.PeriodicSite

        Nsites = len(structure.sites)
        for i in range(Nsites):
            print('\n\n')
            # print(structure.sites[i])
            print(structure.sites[i].species)
            print(structure.sites[i].coords)
            print(structure.sites[i].frac_coords)     
    
    def load_xyz(self, api_key, fractional=False):
        structure = Material().load_structure(api_key, conventional=True)
        Nsites = len(structure.sites)

        xyz_array = np.zeros((Nsites,4))

        if not fractional:
            for i in range(Nsites):
                atom_num = element(str(structure.sites[i].specie)).atomic_number # specie is not a typo
                xyz_array[i] = [atom_num,*structure.sites[i].coords]
        else:
            for i in range(Nsites):
                atom_num = element(str(structure.sites[i].specie)).atomic_number 
                xyz_array[i] = [atom_num,*structure.sites[i].frac_coords]

        return xyz_array

    def to_xyz(self,fractional=False):

        return Material().load_xyz(self.API_KEY, fractional)

    def to_box(self, fractional=False):

        xyz_array = Material().load_xyz(self.API_KEY, fractional)

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

        xyz_array = Material().load_xyz(self.API_KEY, fractional)

        xyz_array[:,1:]*=spacing

        ax = plt.axes(projection='3d')

        colors = np.linspace(0,2**24,118,dtype='int') #divide color range into 118 colors (for the 118 chemical elements)

        for i in xyz_array:
            atom,*xyz = i.astype('int')
            ax.scatter3D(*xyz, s=100, color="#"+hex(colors[atom])[2:])

        set_axes_equal(ax)           

        plt.axis('off')
        plt.show()


    def XRD(self):

        structure = Material().load_structure(self.API_KEY)
        
        # this example shows how to obtain an XRD diffraction pattern
        # these patterns are calculated on-the-fly from the structure
        calculator = XRDCalculator(wavelength='CuKa')
        pattern = calculator.get_pattern(structure)

        print('\npattern:\n',pattern)
        

    def thermo(self):   

        with MPRester(api_key=self.API_KEY) as mpr:

            # for a single material
            # thermo_doc = mpr.thermo.get_data_by_id('mp-1103503')      #   DOES NOT WORK (WHY?)

            # for many materials, it's much faster to use
            # the `search` method, where additional material_ids can 
            # be added to this list
            thermo_docs = mpr.thermo.search(material_ids=[self.structure_ID])
            
        print(thermo_docs)
