from SALib.sample import fast_sampler
import scipy.io as sio
import numpy as np
import math

from sys import argv

def generate_samples(output_dir):

    # Define the model inputs
    problem = {
        'num_vars': 4,
        'names': ['R_filo', 'rho', 'c_i', 'dummy'],
        'bounds': [[30, 75],
                   [0.25, 1],
                   [0, 5],
                   [0, 1]]
    }

    # Generate samples
    #   FAST Sampling with Re-sampling:
    num_samples = 131
    num_resamples = 3
    param_values = fast_sampler.sample(problem, num_samples, seed = 0)
    for i in range(1,num_resamples):
        new_vals = fast_sampler.sample(problem, num_samples, seed = i)
        param_values = np.dstack((param_values, new_vals))

    sio.savemat((output_dir+'eFAST_sample_matrix.mat'), {"eFAST_sample_matrix": param_values})
    print(param_values.shape)

if __name__ == "__main__":
    if len(argv) != 2:
        print("Usage: {} OUTPUT_DIR".format(argv[0]))
        exit(2)

    output_dir = "./"

    generate_samples(output_dir)
