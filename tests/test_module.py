import mlbands

# replace below line with API key [ in string form ] from Materials Project site (https://materialsproject.org/api#api-key)
SECRET_KEY = ''    

def test_materialprops():
    material = mlbands.Material(SECRET_KEY, 'mp-1103503')

    material.structural()
    material.XRD()
    material.thermo()
    totalmag = material.magnetism().total_magnetization
    print(totalmag)

def test_visuals():
    material = mlbands.Material(SECRET_KEY, 'mp-1103502')
    box = material.to_box(True)
    material.visual(10,True)
    box = material.to_box()
    material.visual()

    material = mlbands.Material(SECRET_KEY, 'mp-1103506')
    box = material.to_box(True)
    material.visual(10,True)

def test_loadvisual():
    xdata = mlbands.load('materials.data')
    print(xdata)

    mlbands.Material(SECRET_KEY, box_array = xdata[3]).visual()
    mlbands.Material(SECRET_KEY).visual()


def test_bands():

    training = mlbands.Group(SECRET_KEY)
    training.data_make(range(1,100))
    # training.data_make(range(1,30),True)
    training.resize_boxes()

    testing = mlbands.Group(SECRET_KEY)
    testing.data_make(range(300,350))
    # testing.data_make(range(300,314),True)
    testing.resize_boxes()
    
    machine = mlbands.Machine()
    machine.learn([training.X,training.Y],[testing.X,testing.Y])

    mlbands.save(training, 'train.data')
    mlbands.save(testing, 'test.data')

def test_bands_load():

    training = mlbands.load('train.data')
    trainXY = [training.X,training.Y]
    testing = mlbands.load('test.data')
    testXY = [testing.X,testing.Y]

    machine = mlbands.Machine()
    machine.learn(trainXY,testXY)

def test_dataexpand():

    training = mlbands.Group(SECRET_KEY)
    traindata = mlbands.load('train.data')
    training.transfer(traindata)

    testing = mlbands.load('test.data')
    testdata = mlbands.load('train.data')
    testing.transfer(testdata)

    print(training.X)
    print(testing.materials)


    
test_materialprops()
test_visuals()
test_loadvisual()
test_dataexpand()
# test_bands()
test_bands_load()
