import ase
from ase.io import read,write,Trajectory 
from ase.calculators.vasp import Vasp 
from ase import Atoms
from ase.build import bulk, add_adsorbate
from ase.constraints import FixAtoms
from ase.visualize import view

s = bulk('Si')
# c = FixAtoms(indices = [atom.index for atom in s if atom.position[2] < 14]) # Fixing the layers for which the z coordinate is lesser than 14 (**More the z, upper is the layer)
write('POSCAR',s)
atoms = read('POSCAR', index = -1) 
atoms.pbc = True 
calc=Vasp( 
gga='PE',  # Change this according to the functional used
lreal='Auto',
lplane = True,
lwave = False,
lcharg = False, 
npar = 4, # Change this according to the number of cores used
isym = 2, 
prec='Accurate',
encut = 400, # Change this according to the system
ediff = 1e-6, # Change this according to the speed of simulation required
nelm = 500, 
nelmin = 6,
nelmdl = -9,
algo = 'Very_Fast', 
ismear = -5, # Tetrahedron method
# sigma = 0.05, # For metals, it's 0.2 and for others, it's 0.05    # Sigma doesn't make sense for tetrahedron method.
ibrion = 2, # Conjugate gradient method
ediffg = -0.01, # Change this according to the speed of simulation required
nsw = 500,
potim = 0.15,
isif = 3,
kpts=[13,13,13], # Change this according to the system
gamma=True,   
) 
atoms.set_calculator(calc) 
atoms.get_potential_energy()