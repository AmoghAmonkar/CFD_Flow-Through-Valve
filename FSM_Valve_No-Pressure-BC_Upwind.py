# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 19:08:19 2024

@author: aaamonka
"""
import numpy as np
from tqdm import tqdm
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt 

if __name__ == "__main__":

    #Setting length of the pipe
    #Discretization parameters in x direction
    x_start = 0
    x_end = 12
    dx = 0.0125

    #Creating grids for p, u and v in x dir
    px_grid = np.arange(x_start-(dx/2),x_end+(dx),dx)
    uFx_grid = np.arange(x_start,x_end+(dx),dx)
    vFx_grid = np.arange(x_start-(dx/2),x_end+(dx),dx)

    #Setting diameter of the pipe
    #Discretization parameters in y direction
    y_start = 0
    y_end = 0.6
    dy = 0.0125
    
    #Creating grids for p, u and v in y dir
    py_grid = np.arange(y_start-(dy/2),y_end+(dy),dy)
    uFy_grid = np.arange(y_start-(dy/2),y_end+(dy),dy)
    vFy_grid = np.arange(y_start,y_end+(dy),dy)
    
    #Discretization parameters for t
    t_start = 0
    t_end = 0.002
    dt = 0.00025
    
    #Creating grid for time
    t_grid = np.arange(t_start,t_end,dt)

    #Valve Thickness and Opening
    Valve_Thickness = 0.075
    Nodes_Half_Thick = int(Valve_Thickness/dx)//2
    Valve_Opening = 0.65

    Valve_x_Start = int((len(px_grid))/2) - Nodes_Half_Thick
    Valve_x_End = int((len(px_grid))/2) + Nodes_Half_Thick + 1

    Valve_y_Start = int((len(vFy_grid)-1)*Valve_Opening) 
    Valve_y_End = int((len(vFy_grid))-1) + 1

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
    #Pressure at Inlet
    p_old[0,1:-1] = p_old[1,1:-1]
    #Pressure at Outlet
    p_old[-1,1:-1] = -p_old[-2,1:-1]
    #Bottom Wall
    p_old[1:-1,0] = p_old[1:-1,1]
    #Upper Wall
    p_old[1:-1,-1] = p_old[1:-1,-2]

    #Inlet velocity
    u_in = 1.8
    #Setting inlet velocity
    u_old[0,:] = u_in
    #Top Boundary
    u_old[:,-1] = -u_old[:,-2]
    #Bottom Boundary
    u_old[:,0] = -u_old[:,1]
    #Outlet velocity
    u_old[-1,1:-1] = u_old[-2,1:-1]
    # Apply valve boundary conditions nodes inside valve
    u_old[Valve_x_Start:Valve_x_End, Valve_y_Start+1:Valve_y_End] = 0
    # Apply wall boundary condition on lower edge of the valve
    u_old[Valve_x_Start:Valve_x_End, Valve_y_Start+1] = -u_old[Valve_x_Start:Valve_x_End, Valve_y_Start]

    #Setting inlet velocity
    v_old[0,:] = 0
    #Top Boundary
    v_old[:,-1] = 0
    #Bottom Boundary
    v_old[:,0] = 0
    #Outlet velocity
    v_old[-1,1:-1] = v_old[-2,1:-1]
    # Apply valve boundary conditions nodes inside valve
    v_old[Valve_x_Start+1:Valve_x_End, Valve_y_Start:Valve_y_End] = 0
    # Left Valve Wall boundary conditions
    v_old[Valve_x_Start+1, Valve_y_Start:Valve_y_End] = -v_old[Valve_x_Start, Valve_y_Start:Valve_y_End]
    # Right Valve Wall boundary conditions
    v_old[Valve_x_End-1, Valve_y_Start:Valve_y_End] = -v_old[Valve_x_End, Valve_y_Start:Valve_y_End]

    #Setting Pressure inside valve wall zero
    p_old[Valve_x_Start+1:Valve_x_End, Valve_y_Start+1:Valve_y_End] = 0
    #Setting Pressure Gradient on Bottom Edge Zero
    p_old[Valve_x_Start:Valve_x_End, Valve_y_Start+1] = p_old[Valve_x_Start:Valve_x_End, Valve_y_Start]
    #Setting Pressure Gradient on Left Edge Zero
    p_old[Valve_x_Start+1, Valve_y_Start+1:Valve_y_End] = p_old[Valve_x_Start, Valve_y_Start+1:Valve_y_End]
    #Setting Pressure Gradient on Right Edge Zero
    p_old[Valve_x_End-1, Valve_y_Start+1:Valve_y_End] = p_old[Valve_x_End, Valve_y_Start+1:Valve_y_End]
    #Bottom Wall
    p_old[:,0] = p_old[:,1]
    #Upper Wall
    p_old[:,-1] = p_old[:,-2]
    
    for iter in tqdm(range(len(t_grid))):
        
        CFL = max(u_in, np.max(u_old)) * dt / dx
        if CFL > 1:
            print("CFL condition not satisfied.")
            break
            #dt = dt/CFL
        
        #Calculating u
        #Diffusion in x direction
        diff_x = (nu*((u_old[2:,1:-1] + u_old[:-2,1:-1] - 2*u_old[1:-1,1:-1])/(dx**2) + 
                      (u_old[1:-1,2:] + u_old[1:-1,:-2] - 2*u_old[1:-1,1:-1])/(dy**2)))
        
        #Advection in x direction
        Nu_uu = np.abs(u_old[1:-1,1:-1])/dx
        
        duu_dx = (((u_old[1:-1,1:-1] + u_old[2:,1:-1])**2 - (u_old[1:-1,1:-1] + u_old[:-2,1:-1])**2)/(4*dx) +
                    Nu_uu*(np.abs(u_old[1:-1,1:-1] + u_old[2:,1:-1])*(u_old[1:-1,1:-1] - u_old[1:-1,2:]) - 
                           np.abs(u_old[:-2,1:-1] + u_old[1:-1,1:-1])*(u_old[:-2,1:-1] - u_old[1:-1,1:-1]))/(4*dx))
        
        Nu_vu = np.maximum(np.abs(u_old[1:-1,1:-1])/dx, np.abs(v_old[2:-1,1:])/dy)
        
        dvu_dy = (((v_old[1:-2,1:] + v_old[2:-1,1:])*(u_old[1:-1,1:-1] + u_old[1:-1,2:]) - 
                   (v_old[1:-2,:-1] + v_old[2:-1,:-1])*(u_old[1:-1,:-2] + u_old[1:-1,1:-1]))/(4*dy) +
                    Nu_vu*(np.abs(v_old[1:-2,1:] + v_old[2:-1,1:])*(u_old[1:-1,1:-1] + u_old[1:-1,2:]) - 
                           np.abs(v_old[1:-2,:-1] + v_old[2:-1,:-1])*(u_old[1:-1,:-2] + u_old[1:-1,1:-1]))/(4*dy))
        
        advec_x = duu_dx + dvu_dy
        
        #Explicit Euler
        uF[1:-1,1:-1] = u_old[1:-1,1:-1] + dt*(diff_x - advec_x)
        
        #Boundary Conditions
        #Inlet velocity
        uF[0,1:-1] = u_in
        #Making velocity inside the valve stem 0
        uF[Valve_x_Start:Valve_x_End,Valve_y_Start+1:Valve_y_End] = 0
        #Setting Velocity on nodes left to left vertical valve wall so that velocity on wall is 0
        uF[Valve_x_Start,Valve_y_Start+1:Valve_y_End] = -uF[(Valve_x_Start-1),Valve_y_Start+1:Valve_y_End]
        #Setting Velocity on nodes right to right vertical valve wall so that velocity on wall is 0
        uF[(Valve_x_End-2),Valve_y_Start+1:Valve_y_End] = -uF[(Valve_x_End-1),Valve_y_Start+1:Valve_y_End]
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
        Nu_vv = np.abs(v_old[1:-1,1:-1])/dy
        
        dvv_dy = (((v_old[1:-1,1:-1] + v_old[1:-1,2:])**2 - (v_old[1:-1,1:-1] + v_old[1:-1,:-2])**2)/(4*dy) +
                Nu_vv*(np.abs(v_old[1:-1,1:-1] + v_old[1:-1,2:])*(v_old[1:-1,1:-1] - v_old[1:-1,2:]) - 
                       np.abs(v_old[1:-1,:-2] + v_old[1:-1,1:-1])*(v_old[1:-1,1:-1] - v_old[1:-1,:-2]))/(4*dy))
        
        Nu_uv = np.maximum(np.abs(u_old[1:-1,1:-1])/dx, np.abs(v_old[2:-1,1:])/dy)
        
        duv_dx = (((u_old[:-1,1:-2] + u_old[:-1,2:-1])*(v_old[2:,1:-1] + v_old[1:-1,1:-1]) - 
                   (u_old[1:,1:-2] + u_old[1:,2:-1])*(v_old[1:-1,1:-1] + v_old[:-2,1:-1]))/(4*dx) +
                    Nu_vv*(np.abs(u_old[:-1,1:-2] + u_old[:-1,2:-1])*(v_old[1:-1,1:-1] - v_old[2:,1:-1]) - 
                           np.abs(u_old[1:,1:-2] + u_old[1:,2:-1])*(v_old[:-2,1:-1] - v_old[1:-1,1:-1]))/(4*dx))
        
        advec_y = dvv_dy + duv_dx
        
        #Explicit Euler
        vF[1:-1,1:-1] = v_old[1:-1,1:-1] + dt*(diff_y - advec_y)
        
        #Boundary Conditions
        #Inlet velocity
        vF[0,1:-1] = -vF[0,1:-1]
        #Making velocity inside the valve stem 0
        vF[Valve_x_Start:(Valve_x_End),Valve_y_Start+1:Valve_y_End] = 0
        #Making velocity of the nodes below the bottom wall such that velocity on the wall is zero
        vF[Valve_x_Start:(Valve_x_End),(Valve_y_Start)]= -vF[Valve_x_Start:(Valve_x_End),(Valve_y_Start-1)]
        #Bottom Boundary
        vF[:,0] = 0
        #Top Boundary
        vF[:,-1] = 0
        #Outlet velocity
        vF[-1,1:-1] = vF[-2,1:-1]

        #Computing Div of velocities
        Div_v = (uF[1:,1:-1] - uF[:-1,1:-1])/(2*dx) + (vF[1:-1,1:] - vF[1:-1,:-1])/(2*dy)

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
        #Making velocity inside the valve stem 0
        u_new[Valve_x_Start:(Valve_x_End-1),Valve_y_Start:Valve_y_End] = 0
        #Setting Velocity on nodes left to left vertical valve wall so that velocity on wall is 0
        u_new[Valve_x_Start,Valve_y_Start:Valve_y_End] = -u_new[(Valve_x_Start-1),Valve_y_Start:Valve_y_End]
        #Setting Velocity on nodes right to right vertical valve wall so that velocity on wall is 0
        u_new[(Valve_x_End-2),Valve_y_Start:Valve_y_End] = -u_new[(Valve_x_End-1),Valve_y_Start:Valve_y_End]
        #Bottom Boundary
        u_new[:,0] = -u_new[:,1]
        #Top Boundary
        u_new[:,-1] = -u_new[:,-2]
        #Outlet velocity
        u_new[-1,:] = u_new[-2,:]

        #Boundary Conditions
        #Inlet velocity
        v_new[0,1:-1] = -v_new[0,1:-1]
        #Making velocity inside the valve stem 0
        v_new[Valve_x_Start:(Valve_x_End),Valve_y_Start:Valve_y_End] = 0
        #Making velocity of the nodes below the bottom wall such that velocity on the wall is zero
        v_new[Valve_x_Start:(Valve_x_End),(Valve_y_Start)]= -v_new[Valve_x_Start:(Valve_x_End),(Valve_y_Start-1)]
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
    fig, axs = plt.subplots(3, 1, figsize=(15, 12))
    
    # Pressure field contour
    c1 = axs[0].contourf(X_p, Y_p, p_new[1:-1, 1:-1].T, levels=40, cmap='viridis')
    axs[0].set_title("Pressure Field (p)")
    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")
    fig.colorbar(c1, ax=axs[0], label='Pressure')

    # u-velocity field contour with streamlines
    c2 = axs[1].contourf(X_u, Y_u, u_new[1:, 1:-1].T, levels=40, cmap='jet')
    axs[1].streamplot(X_u, Y_u, u_new[1:, 1:-1].T, v_new[1:-1, 1:].T, color='k', linewidth=0.5, density=0.85)
    axs[1].set_title("u-Velocity Profile")
    axs[1].set_xlabel("x")
    axs[1].set_ylabel("y")
    fig.colorbar(c2, ax=axs[1], label='u-Velocity')

    # v-velocity field contour with streamlines
    c3 = axs[2].contourf(X_v, Y_v, v_new[1:-1, 1:].T, levels=40, cmap='RdBu')
    # Sampling frequency in x and y direction
    stride_x = 4
    stride_y = 40
    #Plotting vectors
    axs[2].quiver(X_u[::stride_x,::stride_y], Y_u[::stride_x,::stride_y], u_new[1:, 1:-1].T[::stride_x, ::stride_y], v_new[1:-1, 1:].T[::stride_x, ::stride_y], color='k', scale=20, scale_units='xy', angles='xy')
    axs[2].set_title("v-Velocity Profile")
    axs[2].set_xlabel("x")
    axs[2].set_ylabel("y")
    fig.colorbar(c3, ax=axs[2], label='v-Velocity')

    # Add rectangular patch for the valve
    valve_width = Valve_Thickness
    valve_height = y_end - (vFy_grid[Valve_y_Start] - y_start)

    # Adding patches on the plots
    for ax in axs:
        valve_patch = Rectangle(
            (uFx_grid[Valve_x_Start], vFy_grid[Valve_y_Start]),
            valve_width,
            valve_height,
            color='white',
            alpha=0.85,
            label='Valve'
        )
        ax.add_patch(valve_patch)

    # Adjust layout
    plt.tight_layout()
    plt.show()