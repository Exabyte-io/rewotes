import pandas as pd  
from rdkit.Chem import AllChem as Chem
import sys
from ast import literal_eval

def find_precision(theoretical, experimental):
	return abs((theoretical-experimental)/experimental)

def find_optimal_basis_sets(smile, precisionInPercent, homoLumoConvergance):
	bestBasisSets = []
	bestHomoLumo = []
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

	if(find_precision(ccpVTZb3lypSum,groundStateSum) <= precisionDecimal):
		bestBasisSets.append("cc-pVTZ + B3LYP")
	if(find_precision(sixthreeoneb3lypsum,groundStateSum) <= precisionDecimal):
		bestBasisSets.append("6-31G + B3LYP")
	if(find_precision(augccpDZb3lypSum,groundStateSum) <= precisionDecimal):
		bestBasisSets.append("aug-cc-pVDZ + B3LYP")
	if((find_precision(ccpVTZb3lypSum,groundStateSum) > precisionDecimal) and (find_precision(sixthreeoneb3lypsum,groundStateSum) > precisionDecimal) and (find_precision(augccpDZb3lypSum,groundStateSum) > precisionDecimal)):
		bestBasisSets.append("cc-pVTZ + B3LYP,6-31G + B3LYP, and aug-cc-pVDZ + B3LYP are too imprecise. Please try another basis set")

	if homoLumoConvergance == "homoLumoConvergance":
		g3HomoLumoSum= 0
		ccpVTZb3lypHomoLumoSum = 0
		augccpDZb3lypHomoLumoSum = 0
		sixthreeoneb3lypHomoLumosum = 0 
		for atom in arrayOfAtomsInSmile:
			g3HomoLumoSum += groundStateDF["Gaussian3 (G3) HOMO-LUMO gap"][atom]	
			ccpVTZb3lypHomoLumoSum += groundStateDF["cc-pVTZ + B3LYP HOMO-LUMO gap"][atom]
			sixthreeoneb3lypHomoLumosum += groundStateDF["6-31G + B3LYP HOMO-LUMO gap"][atom]
			augccpDZb3lypHomoLumoSum += groundStateDF["aug-cc-pVDZ + B3LYP HOMO-LUMO gap"][atom]
		minConvergance = min([g3HomoLumoSum,ccpVTZb3lypHomoLumoSum,augccpDZb3lypHomoLumoSum,sixthreeoneb3lypHomoLumosum])
		if g3HomoLumoSum == minConvergance:
			bestHomoLumo.append("Gaussian3 (G3)")
		if ccpVTZb3lypHomoLumoSum == minConvergance:
			bestHomoLumo.append("cc-pVTZ + B3LYP")
		if sixthreeoneb3lypHomoLumosum == minConvergance:
			bestHomoLumo.append("6-31G + B3LYP")
		if augccpDZb3lypHomoLumoSum == minConvergance:
			bestHomoLumo.append("aug-cc-pVDZ + B3LYP")
	outputJSON = {"Best Basis Sets Based on Ground State Precision": bestBasisSets, "Best Basis Sets Based on HOMO-LUMO Gap Convergance": bestHomoLumo}
	return outputJSON

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print("Syntax: python basissetselector.py [.pdb file name location or smiles string] [percetage of precision wihtout %] [homoLumoConverganceif you want check for the HOMO-LUMO convergance]")
	else:
		print(find_optimal_basis_sets(sys.argv[1], literal_eval(sys.argv[2]), sys.argv[3]))