# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 18:32:24 2024

@author: Amogh Amonkar
"""
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt 

if __name__ == "__main__":

    #Setting length of the pipe
    #Discretization parameters in x direction
    x_start = 0
    x_end = 15
    dx = 0.1
    
    #Creating grids for p, u and v in x dir
    px_grid = np.arange(x_start-(dx/2),x_end+(dx),dx)
    uFx_grid = np.arange(x_start,x_end+(dx),dx)
    vFx_grid = np.arange(x_start-(dx/2),x_end+(dx),dx)

    #Setting diameter of the pipe
    #Discretization parameters in y direction
    y_start = 0
    y_end = 1
    dy = 0.1
    
    #Creating grids for p, u and v in y dir
    py_grid = np.arange(y_start-(dy/2),y_end+(dy),dy)
    uFy_grid = np.arange(y_start-(dy/2),y_end+(dy),dy)
    vFy_grid = np.arange(y_start,y_end+(dy),dy)
    
    #Discretization parameters for t
    t_start = 0
    t_end = 30
    dt = 0.01
    
    #Creating grid for time
    t_grid = np.arange(t_start,t_end,dt)

    #Kinematic Viscosity
    nu = 0.01

    #Inital velocities
    p_old = np.zeros((len(px_grid),len(py_grid)))
    u_old = np.zeros((len(uFx_grid),len(uFy_grid)))
    v_old = np.zeros((len(vFx_grid),len(vFy_grid)))
    
    #Intermediate velocities
    uF = np.zeros_like(u_old)
    vF = np.zeros_like(v_old)
    
    #Corrected velocities
    p_new = np.zeros((len(px_grid),len(py_grid)))
    u_new = np.zeros((len(uFx_grid),len(uFy_grid)))
    v_new = np.zeros((len(vFx_grid),len(vFy_grid)))

    #Pressue Poisson Eqn
    Max_iter = 150
    Convergence_Tol = 0.001

    #Applying boundary conditions on pressure field
    #p_in = 5*98066.5
    #p_old[1,1:-1] = p_in 
    #Pressure at Inlet
    p_old[0,1:-1] = p_old[1,1:-1]
    #Pressure at Outlet
    p_old[-1,1:-1] = -p_old[-2,1:-1]
    #Bottom Wall
    p_old[:,0] = p_old[:,1]
    #Upper Wall
    p_old[:,-1] = p_old[:,-2]

    #Inlet velocity
    u_in = 1
    #Setting inlet velocity
    u_old[0,:] = u_in
    #Bottom Boundary
    u_old[:,0] = -u_old[:,1]
    #Top Boundary
    u_old[:,-1] = -u_old[:,-2]
    #Outlet velocity
    u_old[-1,1:-1] = u_old[-2,1:-1]
    
    for iter in tqdm(range(len(t_grid))):
    #for iter in tqdm(range(12)):
        
        Co_num = max(u_in, np.max(u_old)) * dt / dx
        if Co_num > 1:
            raise ValueError("Courant Number greater than 1")
        
        #Calculating u
        #Diffusion in x direction
        diff_x = (nu*((u_old[2:,1:-1] + u_old[:-2,1:-1] - 2*u_old[1:-1,1:-1])/(dx**2) + 
                      (u_old[1:-1,2:] + u_old[1:-1,:-2] - 2*u_old[1:-1,1:-1])/(dy**2)))
        
        #Advection in x direction
        advec_x = (((u_old[1:-1,1:-1] + u_old[2:,1:-1])**2 - (u_old[1:-1,1:-1] + u_old[:-2,1:-1])**2)/(4*dx) +
                    ((v_old[1:-2,1:] + v_old[2:-1,1:])*(u_old[1:-1,1:-1] + u_old[1:-1,2:]) - 
                     (v_old[1:-2,:-1] + v_old[2:-1,:-1])*(u_old[1:-1,:-2] + u_old[1:-1,1:-1]))/(4*dy))
        
        #Explicit Euler
        uF[1:-1,1:-1] = u_old[1:-1,1:-1] + dt*(diff_x - advec_x)
        
        #Boundary Conditions
        #Inlet velocity
        uF[0,1:-1] = u_in
        #Bottom Boundary
        uF[:,0] = -uF[:,1]
        #Top Boundary
        uF[:,-1] = -uF[:,-2]
        #Outlet velocity
        uF[-1,1:-1] = uF[-2,1:-1]
        
        #Calculating v
        #Diffusion in y direction
        diff_y = (nu*((v_old[2:,1:-1] + v_old[:-2,1:-1] - 2*v_old[1:-1,1:-1])/(dx**2) + 
                      (v_old[1:-1,2:] + v_old[1:-1,:-2] - 2*v_old[1:-1,1:-1])/(dy**2)))
        
        #Advection in y direction
        advec_y = (((v_old[1:-1,1:-1] + v_old[1:-1,2:])**2 - (v_old[1:-1,1:-1] + v_old[1:-1,:-2])**2)/(4*dy) +
                   ((u_old[:-1,1:-2] + u_old[:-1,2:-1])*(v_old[2:,1:-1] + v_old[1:-1,1:-1]) - 
                    (u_old[1:,1:-2] + u_old[1:,2:-1])*(v_old[1:-1,1:-1] + v_old[:-2,1:-1]))/(4*dx))
        
        #Explicit Euler
        vF[1:-1,1:-1] = v_old[1:-1,1:-1] + dt*(diff_y - advec_y)
        
        #Boundary Conditions
        #Inlet velocity
        vF[0,1:-1] = -vF[1,1:-1]
        #Bottom Boundary
        vF[:,0] = 0
        #Top Boundary
        vF[:,-1] = 0
        #Outlet velocity
        vF[-1,1:-1] = vF[-2,1:-1]
        
        #Computing Div of velocities
        Div_v = (uF[1:,1:-1] - uF[:-1,1:-1])/(dx) + (vF[1:-1,1:] - vF[1:-1,:-1])/(dy)

        #b matrix for pressure poisson
        b = Div_v/dt
        
        #Solving pressure poisson equation
        loop_count = 0
        error = 10
        
        while (error >= Convergence_Tol and loop_count < Max_iter):
            
            p_new[1:-1, 1:-1] = (p_old[2:,1:-1] + p_old[:-2,1:-1] + p_old[1:-1,2:] + p_old[1:-1,:-2] - ((dx**2)*b))*0.25
            
            error = np.max(abs(p_new - p_old))
            
            p_old = np.copy(p_new)
            loop_count += 1
        
        #Applying boundary conditions on pressure field
        #Pressure at Inlet
        #p_new[1,1:-1] = p_in 
        p_new[0,1:-1] = p_new[1,1:-1]
        #Bottom Wall
        p_new[:,0] = p_new[:,1]
        #Upper Wall
        p_new[:,-1] = p_new[:,-2]
        #Pressure at Outlet
        p_new[-1,:] = -p_new[-2,:]
        #p_old after enforcing BCs
        p_old = np.copy(p_new)
        
        #Calculating corrected velocities
        u_new[1:,1:-1] = uF[1:,1:-1] - dt*((p_new[2:,1:-1] - p_new[:-2,1:-1]) / (2 * dx))
        v_new[1:-1,1:] = vF[1:-1,1:] - dt*((p_new[1:-1,2:] - p_new[1:-1,:-2]) / (2 * dy))
        
        #Boundary Conditions
        #Inlet velocity
        u_new[0,1:-1] = u_in
        #Bottom Boundary
        u_new[:,0] = -u_new[:,1]
        #Top Boundary
        u_new[:,-1] = -u_new[:,-2]
        #Outlet velocity
        u_new[-1,:] = u_new[-2,:]
        
        #Boundary Conditions
        #Inlet velocity
        v_new[0,1:-1] = -v_new[0,1:-1]
        #Bottom Boundary
        v_new[:,0] = 0
        #Top Boundary
        v_new[:,-1] = 0
        #Outlet velocity
        v_new[-1,1:-1] = v_new[-2,1:-1]
        
        #Replacing old velocities by new for next iteration
        v_old = np.copy(v_new)
        u_old = np.copy(u_new)
    
    # Meshgrid for plotting
    X_p, Y_p = np.meshgrid(px_grid[1:-1], py_grid[1:-1])  # For pressure
    X_u, Y_u = np.meshgrid(uFx_grid[1:], uFy_grid[1:-1])  # For u-velocity
    X_v, Y_v = np.meshgrid(vFx_grid[1:-1], vFy_grid[1:])  # For v-velocity

    # Create subplots
    fig, axs = plt.subplots(3, 1, figsize=(15, 15))
    
    # Pressure field contour
    c1 = axs[0].contourf(X_p, Y_p, p_new[1:-1, 1:-1].T, levels=50, cmap='viridis')
    axs[0].set_title("Pressure Field (p)")
    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")
    fig.colorbar(c1, ax=axs[0], label='Pressure')
    
    # u-velocity field contour with streamlines
    c2 = axs[1].contourf(X_u, Y_u, u_new[1:, 1:-1].T, levels=50, cmap='jet')
    axs[1].streamplot(X_u, Y_u, u_new[1:, 1:-1].T, v_new[1:-1, 1:].T, color='k', linewidth=0.5)
    axs[1].set_title("u-Velocity Profile")
    axs[1].set_xlabel("x")
    axs[1].set_ylabel("y")
    fig.colorbar(c2, ax=axs[1], label='u-Velocity')
    
    # v-velocity field contour with streamlines
    c3 = axs[2].contourf(X_v, Y_v, v_new[1:-1, 1:].T, levels=50, cmap='plasma')
    axs[2].streamplot(X_v, Y_v, u_new[1:, 1:-1].T, v_new[1:-1, 1:].T, color='k', linewidth=0.5)
    axs[2].set_title("v-Velocity Profile")
    axs[2].set_xlabel("x")
    axs[2].set_ylabel("y")
    fig.colorbar(c3, ax=axs[2], label='v-Velocity')

    # Adjust layout
    plt.tight_layout()
    plt.show()