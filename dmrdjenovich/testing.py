from computing import resources
from energy_convergence import EnergyConvergence

dir = "/Users/david/Documents/Internship/Mat3ra/rewotes/dmrdjenovich/example/static-beta-tin"
encut = 11.1
thresh = 0.001
rsx = resources.Resources(1, -1)

conv = EnergyConvergence(dir, encut, thresh, rsx)
conv.run(rsx)
