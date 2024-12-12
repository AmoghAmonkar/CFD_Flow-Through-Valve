# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 11:49:51 2024

@author: Amogh Amonkar
"""
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

# Domain and grid setup
x_start = 0
x_end = 12
dx = 0.5
px_grid = np.arange(x_start - (dx / 2), x_end + dx, dx)
ux_grid = np.arange(x_start, x_end + dx, dx)
vx_grid = np.arange(x_start - (dx / 2), x_end + dx, dx)

y_start = 0
y_end = 0.75
dy = 0.1
py_grid = np.arange(y_start - (dy / 2), y_end + dy, dy)
uy_grid = np.arange(y_start - (dy / 2), y_end + dy, dy)
vy_grid = np.arange(y_start, y_end + dy, dy)

# Valve properties
Valve_Thickness = 0.075
Nodes_Half_Thick = int(Valve_Thickness / dx) // 2
Valve_Opening = 0.65

Valve_x_Start = int((len(px_grid)) / 2) - Nodes_Half_Thick
Valve_x_End = int((len(px_grid)) / 2) + Nodes_Half_Thick + 1

Valve_y_Start = int((len(vy_grid) - 1) * Valve_Opening)
Valve_y_End = int((len(vy_grid)) - 1) + 1

# Domain boundary
x_domain = np.array([x_start, x_end, x_end, x_start, x_start])
y_domain = np.array([y_start, y_start, y_end, y_end, y_start])

PX, PY = np.meshgrid(px_grid, py_grid)
UX, UY = np.meshgrid(ux_grid, uy_grid)
VX, VY = np.meshgrid(vx_grid, vy_grid)

# Plot the domain boundary, point grids, and valve patch
fig, ax = plt.subplots()
ax.plot(x_domain, y_domain, 'k--', label='Domain Boundary', linewidth=2)

ax.plot(PX.flatten(), PY.flatten(), 'g+', label='P Grid Points')
ax.plot(UX.flatten(), UY.flatten(), 'bx', label='U Grid Points')
ax.plot(VX.flatten(), VY.flatten(), 'ro', label='V Grid Points')

# Add rectangular patch for the valve
valve_width = Valve_Thickness
valve_height = y_end - (vy_grid[Valve_y_Start] - y_start)
valve_patch = Rectangle(
    (ux_grid[Valve_x_Start], vy_grid[Valve_y_Start]),
    valve_width,
    valve_height,
    color='gray',
    alpha=0.5,
    label='Valve'
)
ax.add_patch(valve_patch)

# Add legends
ax.legend()

# Labels and show plot
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.show()