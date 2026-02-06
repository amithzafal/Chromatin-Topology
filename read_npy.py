import numpy as np

for npyFile in ["ExtSing_data_matrix.npy","Initial_NoneSing_data_matrix.npy","NoneSing_data_matrix.npy","TopoSing_data_matrix.npy"]:
    print("Input %s" % npyFile)
    
    data = np.load(npyFile, allow_pickle=True)

    print(len(data))

    knot = {}
    for knotState in data:
        #print(knotState)
        knot[knotState[3]] = 0
    for knotState in data:
        #print(knotState)        
        knot[knotState[3]] += 1
    print(knot)
