import pandas as pd  

def find_precision(theoretical, experimental):
	return abs((theoretical-experimental)/experimental)

def find_optimal_basis_sets(atoms, precisionInPercent):
	bestBasisSets = []
	precisionDecimal = precisionInPercent / 100
	groundStateDF = pd.read_csv("groundstateenergies.csv",index_col=0)
	groundStateEnergyForAtom = groundStateDF["Experimental Ground State Energy from AE17 UMN (Hartree)"][atoms]
	if(find_precision(groundStateDF["cc-pVTZ + B3LYP"][atoms],groundStateEnergyForAtom) <= precisionDecimal):
		bestBasisSets.append("cc-pVTZ + B3LYP")
	if(find_precision(groundStateDF["6-31G + B3LYP"][atoms],groundStateEnergyForAtom) <= precisionDecimal):
		bestBasisSets.append("6-31G + B3LYP")
	if(find_precision(groundStateDF["aug-cc-pVDZ + B3LYP"][atoms],groundStateEnergyForAtom) <= precisionDecimal):
		bestBasisSets.append("aug-cc-pVDZ + B3LYP")
	if((find_precision(groundStateDF["aug-cc-pVDZ + B3LYP"][atoms],groundStateEnergyForAtom) > precisionDecimal) and (find_precision(groundStateDF["cc-pVTZ + B3LYP"][atoms],groundStateEnergyForAtom) > precisionDecimal) and (find_precision(groundStateDF["6-31G + B3LYP"][atoms],groundStateEnergyForAtom) > precisionDecimal)):
		bestBasisSets.append("cc-pVTZ + B3LYP,6-31G + B3LYP, and aug-cc-pVDZ + B3LYP are too imprecise. Please try another basis set")

	return bestBasisSets

print(find_optimal_basis_sets("Li", 0.14))