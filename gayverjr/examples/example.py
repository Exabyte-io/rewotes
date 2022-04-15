from BasisSetSelector import read_json

opt,functional,tolerance = read_json('ex2.json')
result = opt.optimize(functional,tolerance,verbose=True)
if result == {}:
    print("None of the selected basis sets satisfy the chosen tolerance.")
else:
    print("Results:")
    for basis,error in result.items():
        print('{}: {}'.format(basis,error))
        