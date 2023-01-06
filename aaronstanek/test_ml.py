import pkg

a = pkg.TrainingDataManager.load_from_numpy_file("materials.npy")

m = pkg.Model(a)

m.double()

print(m.test_standard_deviation(a))

m.train(a, epochs = 10)

print(m.test_standard_deviation(a))
