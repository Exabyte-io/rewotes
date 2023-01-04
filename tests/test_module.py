import mlbands

def test_general():
    material = mlbands.Material()
    material.API_KEY = mlbands.SECRET_KEY

    material.structural()
    # material.XRD()
    # material.thermo()


def test_bands():

    for i in range(1,50):
        material = mlbands.Material()
        material.API_KEY = mlbands.SECRET_KEY
        material.structure_ID = 'mp-'+str(i)
        material.bands()

def test_xyz():
    material = mlbands.Material()
    material.API_KEY = mlbands.SECRET_KEY
    # a=material.to_xyz()
    # b=material.to_xyz(True)
    # print(a,b)

    box = material.to_box(True)
    material.visual(10,True)


    box = material.to_box()
    material.visual()

# test_general()
# test_bands()
test_xyz()