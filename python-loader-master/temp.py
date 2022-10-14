# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from pyMCDS import pyMCDS
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt

mcds = pyMCDS('initial.xml', 
              '/Users/duncanmartinson/Documents/PhysiCell/output')
test = mcds.data['continuum_variables']['oxygen']
xx, yy = mcds.get_2D_mesh()
z_val = 0.00
plane_VEGF = mcds.get_concentrations('oxygen', z_slice=z_val)

num_levels = 21
min_conc = plane_VEGF.min()
max_conc = plane_VEGF.max()
if (max_conc-min_conc)<10**(-6):
    min_conc = 0
    max_conc = 1
my_levels = np.linspace(min_conc, max_conc, num_levels)

fig = plt.figure()
ax1 = plt.subplot(211)
cs = ax1.contourf(xx, yy, plane_VEGF, levels=my_levels)
cbar1 = fig.colorbar(cs, shrink=0.75)
cbar1.set_label('uM')

ax2 = plt.subplot(212)
ax2.contour(xx, yy, plane_VEGF, color='black', levels = my_levels, linewidths=0.5)
cbar1 = fig.colorbar(cs, shrink=0.75)
cbar1.set_label('uM')

cell_df = mcds.get_cell_df() 
cell_df.head()

print(mcds.data['discrete_cells']['total_volume'])

mcds = pyMCDS('output00000001.xml', 
              '/Users/duncanmartinson/Documents/PhysiCell/output')
test = mcds.data['continuum_variables']['oxygen']
xx, yy = mcds.get_2D_mesh()
z_val = 0.00
plane_VEGF = mcds.get_concentrations('oxygen', z_slice=z_val)

num_levels = 21
min_conc = plane_VEGF.min()
max_conc = plane_VEGF.max()
if (max_conc-min_conc)<10**(-6):
    min_conc = 0
    max_conc = 1
my_levels = np.linspace(min_conc, max_conc, num_levels)

fig3 = plt.figure()
ax1 = plt.subplot(211)
cs = ax1.contourf(xx, yy, plane_VEGF, levels=my_levels)
cbar1 = fig3.colorbar(cs, shrink=0.75)
cbar1.set_label('uM')

ax2 = plt.subplot(212)
ax2.contour(xx, yy, plane_VEGF, color='black', levels = my_levels, linewidths=0.5)
cbar1 = fig3.colorbar(cs, shrink=0.75)
cbar1.set_label('uM')

print(mcds.data['discrete_cells']['total_volume'])
"""
mcds = pyMCDS('output00000002.xml', 
              '/Users/duncanmartinson/Documents/PhysiCell/output')
test = mcds.data['continuum_variables']['oxygen']
xx, yy = mcds.get_2D_mesh()
z_val = 0.00
plane_VEGF = mcds.get_concentrations('oxygen', z_slice=z_val)

num_levels = 21
min_conc = plane_VEGF.min()
max_conc = plane_VEGF.max()
my_levels = np.linspace(min_conc, max_conc, num_levels)

fig2 = plt.figure()
ax1 = plt.subplot(211)
cs = ax1.contourf(xx, yy, plane_VEGF, levels=my_levels)
cbar1 = fig2.colorbar(cs, shrink=0.75)
cbar1.set_label('uM')

ax2 = plt.subplot(212)
ax2.contour(xx, yy, plane_VEGF, color='black', levels = my_levels, linewidths=0.5)
cbar1 = fig2.colorbar(cs, shrink=0.75)
cbar1.set_label('uM')

print(mcds.data['discrete_cells']['total_volume'])

"""
