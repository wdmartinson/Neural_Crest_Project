#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 17:21:49 2020

@author: duncanmartinson
"""

from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt
 
# load data
mcds = pyMCDS('output00003696.xml', '/Users/duncanmartinson/Documents/PhysiCell/output')
 
# Set our z plane and get our substrate values along it
z_val = 0.00
 
# get our cells data and figure out which cells are in the plane
cell_df = mcds.get_cell_df()
ds = mcds.get_mesh_spacing()
inside_plane = (cell_df['position_z'] < z_val + ds) & (cell_df['position_z'] > z_val - ds)
plane_cells = cell_df[inside_plane]
 
# We're going to plot two types of cells and we want it to look nice
colors = ['black', 'grey']
sizes = [20, 8]
labels = ['Alive', 'Dead']
 
# set up the figure area and add microenvironment layer
fig, ax = plt.subplot()
cs = ax.contourf(xx, yy, plane_oxy, levels=my_levels)
 
# get our cells of interest
alive_cells = plane_cells[plane_cells['cycle_model'] < 6]
dead_cells = plane_cells[plane_cells['cycle_model'] > 6]
 
# plot the cell layer
for i, plot_cells in enumerate((alive_cells, dead_cells)):
    ax.scatter(plot_cells['position_x'].values, 
            plot_cells['position_y'].values, 
            facecolor='none', 
            edgecolors=colors[i],
            alpha=0.6,
            s=sizes[i],
            label=labels[i])
 
# Now we need to add our color bar
cbar1 = fig.colorbar(cs, shrink=0.75)
cbar1.set_label('mmHg')
 
# Let's put the time in to make these look nice
ax.set_aspect('equal')
ax.set_xlabel('x (micron)')
ax.set_ylabel('y (micron)')
ax.set_title('oxygen (mmHg) at t = {:.1f} {:s}, z = {:.2f} {:s}'.format(
                                        mcds.get_time(),
                                        mcds.data['metadata']['time_units'],
                                        z_val,
                                        mcds.data['metadata']['spatial_units'])
ax.legend(loc='upper right')
 
plt.show()