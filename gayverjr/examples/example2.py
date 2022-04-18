from BasisSetSelector import BasisSetOptimizer

opt = BasisSetOptimizer(['STO-3G','6-31G','3-21G','6-31+G*','pc-2'],'homo lumo gap')
# B3LYP/aug-cc-pVTZ result
opt.add_molecule('uracil.xyz',3.864)
result = opt.optimize('b3lyp',0.1, verbose=True)
if result == {}:
    print("None of the selected basis sets satisfy the chosen tolerance.")
else:
    print("Results:")
    for basis,error in result.items():
        print('{}: {}'.format(basis,error))
        