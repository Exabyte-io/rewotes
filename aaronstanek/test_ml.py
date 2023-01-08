import pkg

a = pkg.TrainingDataManager.load_from_numpy_file('materials.npy')

m = pkg.Model(a)

m.double()

print(m.test_standard_deviation(a))
print(m.test_with_error(a, iteration_count=50))

m.train(a, epochs=10)

print(m.test_standard_deviation(a))
print(m.test_with_error(a, iteration_count=50))
