from computing import resources
import energy_convergence

dir = "/Users/david/Documents/Internship/Mat3ra/rewotes/dmrdjenovich/example/static-beta-tin"
encut = 11.1
thresh = 0.001
rsx = resources.Resources(1, -1)

conv = energy_convergence.EnergyConvergence(dir, encut, thresh, rsx)

