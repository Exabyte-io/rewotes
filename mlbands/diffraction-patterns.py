from mp_api.client import MPRester
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import secret

with MPRester(api_key=secret.API_KEY) as mpr:
    # first retrieve the relevant structure
    structure = mpr.get_structure_by_material_id('mp-1103503')

# important to use the conventional structure to ensure
# that peaks are labelled with the conventional Miller indices
sga = SpacegroupAnalyzer(structure)
conventional_structure = sga.get_conventional_standard_structure()

# this example shows how to obtain an XRD diffraction pattern
# these patterns are calculated on-the-fly from the structure
calculator = XRDCalculator(wavelength='CuKa')
pattern = calculator.get_pattern(conventional_structure)

print(pattern)