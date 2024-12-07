# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 11:49:51 2024

@author: Amogh Amonkar
"""
import numpy as np
import matplotlib.pyplot as plt

x_start = 0
x_end = 15
dx = 0.5
px_grid = np.arange(x_start-(dx/2),x_end+(dx),dx)
ux_grid = np.arange(x_start,x_end+(dx),dx)
vx_grid = np.arange(x_start-(dx/2),x_end+(dx),dx)

y_start = 0
y_end = 1
dy = 0.5
py_grid = np.arange(y_start-(dy/2),y_end+(dy),dy)
uy_grid = np.arange(y_start-(dy/2),y_end+(dy),dy)
vy_grid = np.arange(y_start,y_end+(dy),dy)

x_domain = np.array([x_start,x_end,x_end,x_start,x_start])
y_domain = np.array([y_start,y_start,y_end,y_end,y_start])

PX, PY = np.meshgrid(px_grid,py_grid)
UX, UY = np.meshgrid(ux_grid,uy_grid)
VX, VY = np.meshgrid(vx_grid,vy_grid)

# Plot the domain boundary, point grids, and vector grids
plt.plot(x_domain, y_domain,'k--', label='Domain Boundary',linewidth=2)
plt.plot(PX.flatten(), PY.flatten(),'g+', label='P Grid Points')
plt.plot(UX.flatten(), UY.flatten(),'bx', label='U Grid Points')
plt.plot(VX.flatten(), VY.flatten(),'ro', label='V Grid Points')

# Add legends
plt.legend()

# Labels and show plot
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.grid()
plt.show()