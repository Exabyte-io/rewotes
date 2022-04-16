from BasisSetSelector import read_json

opt = read_json('ex2.json')
opt.add_basis_set('sto-3g')
opt.remove_basis_set('6-31G')
result = opt.optimize('b3lyp',0.5,verbose=True)
if result == {}:
    print("None of the selected basis sets satisfy the chosen tolerance.")
else:
    print("Results:")
    for basis,error in result.items():
        print('{}: {}'.format(basis,error))
        