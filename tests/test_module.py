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

def test_loadvisual():
    xdata = mlbands.load('materials.data')
    print(xdata)

    mlbands.Material(mlbands.SECRET_KEY, box_array = xdata[3]).visual()
    mlbands.Material(mlbands.SECRET_KEY).visual()


def test_bands():

    training = mlbands.Group(mlbands.SECRET_KEY)
    training.make_data(range(1,30),True)
    training.resize()

    testing = mlbands.Group(mlbands.SECRET_KEY)
    testing.make_data(range(300,314),True)
    testing.resize()
    
    machine = mlbands.Machine()
    machine.learn([training.X,training.Y],[testing.X,testing.Y])

    # mlbands.save([training.X,training.Y], 'train.data')
    # mlbands.save([testing.X,testing.Y], 'test.data')

def test_bands_load():

    train = mlbands.load('train.data')
    test = mlbands.load('test.data')
    machine = mlbands.Machine()
    machine.learn(train,test)



    
# test_material()
# test_visuals()
# test_loadvisual()
test_bands()
# test_bands_load()