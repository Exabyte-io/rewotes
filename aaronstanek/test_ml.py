import pkg

a = pkg.TrainingDataManager.load_from_archive_file("materials.materialarchive")

m = pkg.Model(a)

m.double()

print(m.test_standard_deviation(a))

m.train(a, epochs = 10)

print(m.test_standard_deviation(a))
