import mlbands

def test_general():
    material = mlbands.Material()
    material.API_KEY = mlbands.SECRET_KEY

    material.structural()
    material.XRD()
    material.thermo()


def test_bands():

    for i in range(1,50):
        material = mlbands.Material()
        material.API_KEY = mlbands.SECRET_KEY
        material.structure_ID = 'mp-'+str(i)
        material.bands()

test_general()
test_bands()