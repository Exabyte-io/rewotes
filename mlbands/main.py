from mp_api.client import MPRester
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class Material:
    def __init__(self):
        self.API_KEY = ''
        self.structure_ID = 'mp-1103503'

    def bands(self):
        with MPRester(api_key=self.API_KEY) as mpr:
            bandstructure = mpr.get_bandstructure_by_material_id(self.structure_ID)

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

        # with MPRester(api_key=self.API_KEY) as mpr:
        #     # first retrieve the relevant structure
        #     structure = mpr.get_structure_by_material_id(self.structure_ID)

        conventional_structure = Material().load_structure(self.API_KEY)


        # print('\nstructure:\n{}\n\nsga:\n{}\n\nconventional structure:\n{}'.\
        #     format(structure,sga,conventional_structure))

        # print(conventional_structure)

        print(conventional_structure.lattice)
        print(conventional_structure.sites)  #https://pymatgen.org/pymatgen.core.sites.html?highlight=periodicsite#pymatgen.core.sites.PeriodicSite

        Nsites = len(conventional_structure.sites)
        for i in range(Nsites):
            print('\n\n')
            # print(conventional_structure.sites[i])
            print(conventional_structure.sites[i].species)
            print(conventional_structure.sites[i].coords)
            print(conventional_structure.sites[i].frac_coords)     

    def XRD(self):

        conventional_structure = Material().load_structure(self.API_KEY)
        
        # this example shows how to obtain an XRD diffraction pattern
        # these patterns are calculated on-the-fly from the structure
        calculator = XRDCalculator(wavelength='CuKa')
        pattern = calculator.get_pattern(conventional_structure)

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
