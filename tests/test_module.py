import mlbands

def test_general():
    material = mlbands.Material()
    material.API_KEY = mlbands.SECRET_KEY

    material.structural()
    material.XRD()
    material.thermo()



test_general()
