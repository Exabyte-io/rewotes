import pandas as pd  
from rdkit.Chem import AllChem as Chem
# from rdkit.Chem.rdmolfiles import MolFromPDB


def find_precision(theoretical, experimental):
	return abs((theoretical-experimental)/experimental)

def find_optimal_basis_sets(smile, precisionInPercent,):
	bestBasisSets = []
	precisionDecimal = precisionInPercent / 100
	groundStateDF = pd.read_csv("groundstateenergies.csv",index_col=0)
	if(".pdb" in smile):
		mol = Chem.MolFromPDBFile(smile)
	else:	
		mol = Chem.MolFromSmiles(smile)
	arrayOfAtomsInSmile = [atom.GetSymbol() for atom in mol.GetAtoms()]
	groundStateSum = 0
	ccpVTZb3lypSum = 0
	augccpDZb3lypSum = 0
	sixthreeoneb3lypsum = 0
	for atom in arrayOfAtomsInSmile:
		groundStateSum += groundStateDF["Experimental Ground State Energy from AE17 UMN (Hartree)"][atom]	
		ccpVTZb3lypSum += groundStateDF["cc-pVTZ + B3LYP"][atom]
		sixthreeoneb3lypsum += groundStateDF["6-31G + B3LYP"][atom]
		augccpDZb3lypSum += groundStateDF["aug-cc-pVDZ + B3LYP"][atom]
	print(sixthreeoneb3lypsum)
	if(find_precision(ccpVTZb3lypSum,groundStateSum) <= precisionDecimal):
		bestBasisSets.append("cc-pVTZ + B3LYP")
	if(find_precision(sixthreeoneb3lypsum,groundStateSum) <= precisionDecimal):
		bestBasisSets.append("6-31G + B3LYP")
	if(find_precision(augccpDZb3lypSum,groundStateSum) <= precisionDecimal):
		bestBasisSets.append("aug-cc-pVDZ + B3LYP")
	if((find_precision(ccpVTZb3lypSum,groundStateSum) > precisionDecimal) and (find_precision(sixthreeoneb3lypsum,groundStateSum) > precisionDecimal) and (find_precision(augccpDZb3lypSum,groundStateSum) > precisionDecimal)):
		bestBasisSets.append("cc-pVTZ + B3LYP,6-31G + B3LYP, and aug-cc-pVDZ + B3LYP are too imprecise. Please try another basis set")

	return bestBasisSets

print(find_optimal_basis_sets("proline.pdb", 0.01))