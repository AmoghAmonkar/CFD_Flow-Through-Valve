import numpy as np
import matplotlib.pyplot as plt

x_start = 0
x_end = 15
dx = 0.05
px_grid = np.arange(x_start-(dx/2),x_end+(dx),dx)
ux_grid = np.arange(x_start,x_end+(dx),dx)
vx_grid = np.arange(x_start-(dx/2),x_end+(dx),dx)

y_start = 0
y_end = 1
dy = 0.05
py_grid = np.arange(y_start-(dy/2),y_end+(dy),dy)
uy_grid = np.arange(y_start-(dy/2),y_end+(dy),dy)
vy_grid = np.arange(y_start,y_end+(dy),dy)

#Valve Thickness and Opening
Valve_Thickness = 0.2
Nodes_Half_Thick = int(Valve_Thickness/dx)//2
Valve_Opening = 0.1

Valve_x_Start = int((len(px_grid))/2) - Nodes_Half_Thick
Valve_x_End = int((len(px_grid))/2) + Nodes_Half_Thick + 1

Valve_y_Start = int((len(py_grid)-2)*Valve_Opening) 
Valve_y_End = int((len(py_grid))-2)

x_domain = np.array([x_start,x_end,x_end,x_start,x_start])
y_domain = np.array([y_start,y_start,y_end,y_end,y_start])

#print(px_grid[Valve_x_Start],px_grid[Valve_x_End])
#print(py_grid[Valve_y_Start],py_grid[Valve_y_End])

PX, PY = np.meshgrid(px_grid,py_grid)
#UX, UY = np.meshgrid(ux_grid,uy_grid)
#VX, VY = np.meshgrid(vx_grid,vy_grid)
Valve_x, Valve_y = np.meshgrid(px_grid[Valve_x_Start:Valve_x_End],py_grid[Valve_y_Start:Valve_y_End])

# Plot the domain boundary, point grids, and vector grids
plt.plot(x_domain, y_domain,'k--', label='Domain Boundary',linewidth=2)
#plt.plot(PX.flatten(), PY.flatten(),'g+', label='P Grid Points')
#plt.plot(UX.flatten(), UY.flatten(),'bx', label='U Grid Points')
#plt.plot(VX.flatten(), VY.flatten(),'ro', label='V Grid Points')
plt.plot(Valve_x.flatten(), Valve_y.flatten(),'rx', label='Valve Boundary',linewidth=0.5)

# Add legends
plt.legend()

# Labels and show plot
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.grid()
plt.show()