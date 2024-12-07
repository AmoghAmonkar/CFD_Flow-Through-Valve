# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 20:28:41 2024

@author: aaamonka
"""
import numpy as np

if __name__ == "__main__":

    #Setting length of the pipe
    #Discretization parameters in x direction
    x_start = 0
    x_end = 15
    dx = 0.05

    #Creating grids for p, u and v in x dir
    px_grid = np.arange(x_start-(dx/2),x_end+(dx),dx)
    uFx_grid = np.arange(x_start,x_end+(dx),dx)
    vFx_grid = np.arange(x_start-(dx/2),x_end+(dx),dx)

    #Setting diameter of the pipe
    #Discretization parameters in y direction
    y_start = 0
    y_end = 1
    dy = 0.05
    
    #Creating grids for p, u and v in y dir
    py_grid = np.arange(y_start-(dy/2),y_end+(dy),dy)
    uFy_grid = np.arange(y_start-(dy/2),y_end+(dy),dy)
    vFy_grid = np.arange(y_start,y_end+(dy),dy)
    
    #Discretization parameters for t
    t_start = 0
    t_end = 0.105
    dt = 0.001
    
    #Creating grid for time
    t_grid = np.arange(t_start,t_end,dt)

    #Valve Thickness and Opening
    Valve_Thickness = 0.2
    Nodes_Half_Thick = int(Valve_Thickness/dx)//2
    Valve_Opening = 0.8

    Valve_x_Start = int((len(px_grid))/2) - Nodes_Half_Thick
    Valve_x_End = int((len(px_grid))/2) + Nodes_Half_Thick + 1

    Valve_y_Start = int((len(vFy_grid)-1)*Valve_Opening) 
    Valve_y_End = int((len(vFy_grid))-1) + 1

    #Kinematic Viscosity
    nu = 0.01

    #Inital velocities
    p_old = np.ones((len(px_grid),len(py_grid)))
    u_old = np.ones((len(uFx_grid),len(uFy_grid)))
    v_old = np.ones((len(vFx_grid),len(vFy_grid)))
    
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
    p_old[:,0] = p_old[:,1]
    #Upper Wall
    p_old[:,-1] = p_old[:,-2]

    #Inlet velocity
    u_in = 1.8
    #Setting inlet velocity
    u_old[0,:] = u_in
    # Apply valve boundary conditions nodes inside valve
    u_old[Valve_x_Start:Valve_x_End, Valve_y_Start+1:Valve_y_End] = 0
    # Apply wall boundary condition on lower edge of the valve
    u_old[Valve_x_Start:Valve_x_End, Valve_y_Start+1] = -u_old[Valve_x_Start:Valve_x_End, Valve_y_Start]
    #Top Boundary
    u_old[:,-1] = -u_old[:,-2]
    #Bottom Boundary
    u_old[:,0] = -u_old[:,1]
    #Outlet velocity
    u_old[-1,1:-1] = u_old[-2,1:-1]

    #Setting inlet velocity
    v_old[0,:] = 0
    # Apply valve boundary conditions nodes inside valve
    v_old[Valve_x_Start+1:Valve_x_End, Valve_y_Start:Valve_y_End] = 0
    # Left Valve Wall boundary conditions
    v_old[Valve_x_Start+1, Valve_y_Start:Valve_y_End] = -v_old[Valve_x_Start, Valve_y_Start:Valve_y_End]
    # Right Valve Wall boundary conditions
    v_old[Valve_x_End-1, Valve_y_Start:Valve_y_End] = -v_old[Valve_x_End, Valve_y_Start:Valve_y_End]
    #Top Boundary
    v_old[:,-1] = 0
    #Bottom Boundary
    v_old[:,0] = 0
    #Outlet velocity
    v_old[-1,1:-1] = v_old[-2,1:-1]

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