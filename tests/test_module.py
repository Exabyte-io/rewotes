import mlbands

def test_general():
    material = mlbands.Material()
    material.API_KEY = mlbands.SECRET_KEY

    material.structural()
    material.XRD()
    material.thermo()


def test_bands():
    material = mlbands.Material()
    material.API_KEY = mlbands.SECRET_KEY

    material.bands()

test_general()
test_bands()