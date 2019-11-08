def readmatfile(datapath, filename,listvar):
    print('#############################################################################')
    # datapath = '...' path to folder where mat file is stored
    # filename = 'file.mat' is the name of the dataset
    # listvar = boolean. If = 1, lists out the info about each variable, if = 0 this info is hidden
    import numpy as np
    import sys
    sys.path.insert(0, datapath)
    import hdf5storage
    data = hdf5storage.loadmat(datapath + filename)

##    # An alternative if you have a matlab engine- takes longer to load
##    import matlab.engine
##    eng = matlab.engine.start_matlab()
##    matdata = eng.load(datapath + filename,nargout=1)
##    
    pydata = {}
    if listvar == 1:
        for key in matdata:
            global value
            value = matdata[key]
            value = np.asarray(value)
            pydata[key] = value
            print("Name:",key, "  |  Size:", value.shape, "  |  Type:", type(value))
    print('#############################################################################')
    return pydata

    
