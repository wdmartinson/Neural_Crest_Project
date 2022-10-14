from SALib.analyze import fast
import scipy.io as sio
import numpy as np

problem = {
    'num_vars': 4,
    'names': ['R_filo', 'rho', 'c_i', 'dummy'],
    'bounds': [[30, 75],
               [0.25, 1],
               [0, 5],
               [0, 1]]
}

param_values = sio.loadmat('eFAST_sample_matrix.mat')['eFAST_sample_matrix']
output = sio.loadmat('eFAST_output_stats.mat')

dicts = {}
keys = list(output.keys())[3:]
for i in range(len(keys)):
    name = keys[i]
    Y = output[name]
    Sob_indices = []
    for j in range(param_values.shape[2]):
        # Perform analysis
        Z = Y[:,0,j]
        Si = fast.analyze(problem, Z, print_to_console = True)
        Sob_indices.append(Si)
    dicts[name] = Sob_indices
    sio.savemat(name+'_sobol_indices.mat', {'sobol_indices': Sob_indices})
keys2 = list(dicts.keys())
print("\n")
for ii in range(len(keys2)):
    name = keys2[ii]
    print("Sobol Indices, "+name+":")
    for jj in range(len(dicts[name])):
        print(dicts[name][jj]['S1'])
        print(dicts[name][jj]['ST'])
    print("\n")
