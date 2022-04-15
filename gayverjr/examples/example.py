from BasisSetSelector import read_json

opt,functional,tolerance = read_json('ex1.json')
result = opt.optimize(functional,tolerance)
if result == {}:
    print("None of the selected basis sets satisfy the chosen tolerance.")
else:
    print("Results:")
    for basis,error in result.items():
        print('{}: {}'.format(basis,error))
        