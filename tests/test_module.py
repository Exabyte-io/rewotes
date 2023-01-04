import mlbands

def test_general():
    material = mlbands.Material('mp-1103503')
    material.API_KEY = mlbands.SECRET_KEY

    material.structural()
    material.XRD()
    material.thermo()


def test_bands():

    for i in range(1,10):
        material = mlbands.Material('mp-1103503')
        material.API_KEY = mlbands.SECRET_KEY
        material.material_ID = 'mp-'+str(i)
        material.bands()

def test_xyz():
    material = mlbands.Material('mp-1103502')
    material.API_KEY = mlbands.SECRET_KEY

    box = material.to_box(True)
    material.visual(10,True)

    box = material.to_box()
    material.visual()

# test_general()
# test_bands()
test_xyz()