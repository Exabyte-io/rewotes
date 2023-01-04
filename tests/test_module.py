import mlbands

def test_material():
    material = mlbands.Material(mlbands.SECRET_KEY, 'mp-1103503')

    material.structural()
    material.XRD()
    material.thermo()


def test_visuals():
    material = mlbands.Material(mlbands.SECRET_KEY, 'mp-1103502')
    box = material.to_box(True)
    material.visual(10,True)
    box = material.to_box()
    material.visual()

    material = mlbands.Material(mlbands.SECRET_KEY, 'mp-1103506')
    box = material.to_box(True)
    material.visual(10,True)



def test_bands():

    group = mlbands.Group(mlbands.SECRET_KEY)
    data = group.make_data()

    print(data)
    
test_material()
test_visuals()
test_bands()