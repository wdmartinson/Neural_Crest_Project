# run the following commands in the command line/terminal, if ripser is not
#    already installed:
# pip install cython
# pip install ripser
from ripser import ripser
# from persim import plot_diagrams as plot_dgms
import numpy as np

# Load in the positions of the cells at the final time point:
points = np.loadtxt('final_cell_positions_puncta_model.csv', delimiter = ',')
# Use ripser to compute persistence homology diagrams, with Rips complexes:
diagrams = ripser(points,maxdim=0,thresh=500.0)['dgms']
# plot_dgms(diagrams)

# Save the information into another csv file, which you can write to in matlab
np.savetxt('persistence_homology_diagrams_puncta_model.csv', diagrams[0], delimiter=',')
