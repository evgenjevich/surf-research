
# coding: utf-8

# In[ ]:

import numpy as np
import sympy
import fipy as fp
from fipy import numerix as nmx
import matplotlib.pyplot as plt
import os
import sys

c, rho_s, c_alpha, c_beta = sympy.symbols("c_var rho_s c_alpha c_beta")
f_0 = rho_s * (c - c_alpha)**2 * (c_beta - c)**2

sympy.diff(f_0, c, 2)

# command format:
# python 1_edit.py a 2.0 100
# where nx = 100
# and dx = 2.0
# for benchmark problem 1a
# or
# python 1_edit.py a 100 2500
# for a sphere, radius 100 and containing 2500 cells

problem = sys.argv[1]

if (problem == 'a') or (problem == 'b') or (problem == 'c'):
    nx = int(sys.argv[3])
    dx = float (sys.argv[2])
    epsilon = .01

if (problem == "a"):
    print "Creating mesh for problem a"
    mesh = fp.PeriodicGrid2D(nx=nx, ny=nx, dx=dx, dy=dx)
    
elif (problem == "b"):
    print "Creating mesh for problem b"
    mesh = fp.Grid2D(nx=nx, ny=nx, dx=dx, dy=dx)
    
elif (problem == "c"):
    print "Creating mesh for problem c"
#    mesh = fp.Grid2D(dx=0.5, dy=0.5, nx=40, ny=200) + (fp.Grid2D(dx=0.5, dy=0.5, nx=200, ny=40) + [[-40],[100]])
    mesh = fp.Grid2D(Lx=20., Ly=100.0, nx=nx / 5, ny=nx) + (fp.Grid2D(Ly=20.0, Lx=100.0, nx=nx, ny=nx / 5) + [[-40],[100]])
       
elif (problem == "d"):
    print "Creating mesh for problem d"
    r = float(sys.argv[2])
    numCellsDesired = int(sys.argv[3])
    epsilon = .05
    cellSize = 16 * np.pi * r**2 / (nmx.sqrt(3.) * numCellsDesired)
    cellSize = nmx.sqrt(cellSize)
    
    substring1 = '''
    radius = {0};
    cellSize = {1};
    '''.format(r, round(cellSize, 6))

    mesh = fp.Gmsh2DIn3DSpace(substring1 + '''

    // create inner 1/8 shell
    Point(1) = {0, 0, 0, cellSize};
    Point(2) = {-radius, 0, 0, cellSize};
    Point(3) = {0, radius, 0, cellSize};
    Point(4) = {0, 0, radius, cellSize};
    Circle(1) = {2, 1, 3};
    Circle(2) = {4, 1, 2};
    Circle(3) = {4, 1, 3};
    Line Loop(1) = {1, -3, 2};
    Ruled Surface(1) = {1};

    // create remaining 7/8 inner shells
    t1[] = Rotate {{0,0,1},{0,0,0},Pi/2} {Duplicata{Surface{1};}};
    t2[] = Rotate {{0,0,1},{0,0,0},Pi} {Duplicata{Surface{1};}};
    t3[] = Rotate {{0,0,1},{0,0,0},Pi*3/2} {Duplicata{Surface{1};}};
    t4[] = Rotate {{0,1,0},{0,0,0},-Pi/2} {Duplicata{Surface{1};}};
    t5[] = Rotate {{0,0,1},{0,0,0},Pi/2} {Duplicata{Surface{t4[0]};}};
    t6[] = Rotate {{0,0,1},{0,0,0},Pi} {Duplicata{Surface{t4[0]};}};
    t7[] = Rotate {{0,0,1},{0,0,0},Pi*3/2} {Duplicata{Surface{t4[0]};}};

    // create entire inner and outer shell
    Surface Loop(100)={1, t1[0],t2[0],t3[0],t7[0],t4[0],t5[0],t6[0]};
    ''', order=2.0).extrude(extrudeFunc=lambda r: 1.01*r)

c_alpha = 0.3
c_beta = 0.7
kappa = 2.0
M = 5.0
c_0 = 0.5
rho_s = 5.0

dt_max = 1.0

# solution variable
c_var = fp.CellVariable(mesh=mesh, name=r"$c$", hasOld=True)

# array of sample c-values: used in f versus c plot
vals = np.linspace(-.1, 1.1, 1000)

if (problem == 'a' or 'b' or 'c'):
    # 2D mesh coordinates
    x, y = np.array(mesh.x), np.array(mesh.y)
    # initial value for square and T domains
    c_var[:] = c_0 + epsilon * (np.cos(0.105 * x) * np.cos(0.11 * y) + (np.cos(0.13 * x) * np.cos(0.087 * y))**2 + np.cos(0.025 * x - 0.15 * y) * np.cos(0.07 * x - 0.02 * y))
if (problem == 'd'):
    print "number of cells: " , mesh.numberOfCells
    # 3D mesh coordinates
    x, y, z = np.array(mesh.x), np.array(mesh.y), np.array(mesh.z)
    
    # convert from rectangular to spherical coordinates
    theta = fp.CellVariable(name=r"$\theta$", mesh=mesh)
    theta = nmx.arctan2(z, nmx.sqrt(x**2 + y**2))
    phi = fp.CellVariable(name=r"$\phi$", mesh=mesh)
    phi = nmx.arctan2(y, x)
     
    # initial value for spherical domain
    c_var[:]  = c_0 + epsilon * ((np.cos(8*theta))*(np.cos(15*phi)) + ((np.cos(12*theta))*(np.cos(10*phi)))**2 + ((np.cos(2.5*theta - 1.5*phi))*(np.cos(7*theta - 2*phi))))

# bulk free energy density
def f_0(c):
    return rho_s*((c - c_alpha)**2)*((c_beta-c)**2)
def f_0_var(c_var):
    return 2*rho_s*((c_alpha - c_var)**2 + 4*(c_alpha - c_var)*(c_beta - c_var) + (c_beta - c_var)**2)
# free energy
def f(c):
    return (f_0(c)+ .5*kappa*(c.grad.mag)**2)

# plot free energy density versus c
def plotf_c():
    plt.figure(1)
    plt.xlabel('c')
    plt.ylabel('f_0')
    plt.plot(vals, f_0(vals))
    plt.show()
    
# save elapsed time and free energy at each data point
time_data = []
cvar_data = []
f_data = []
# checks whether a folder for the data from this simulation exists
# if not, creates one in the home directory
file_dir = "/data/and9/surf-research/Anushka/1{0}".format(problem)
if not os.path.exists(file_dir):
    os.makedirs(file_dir)
def save_data(time, cvar, f, step):
    time_data.append(time)
    cvar_data.append(np.array(cvar.value))
    f_data.append(f.value)
    
    file_name = file_dir + "/1{0}{1}_{2}".format(problem, int(sys.argv[3]), str(step).rjust(5, '0'))
    np.savez(file_name, time = time_data, c_var = cvar_data, f = f_data)
    

# solver equation    
eqn = fp.TransientTerm(coeff=1.) == fp.DiffusionTerm(M * f_0_var(c_var)) - fp.DiffusionTerm((M, kappa))

elapsed = 0.0
steps = 0
dt = 0.01
total_sweeps = 2
tolerance = 1e-1

# controls on how long the simulation runs: steps, duration, or both
total_steps = 20
duration = 10000

c_var.updateOld()
from fipy.solvers.pysparse import LinearLUSolver as Solver
solver = Solver()
print "Starting Solver."
while elapsed <= duration:
    res0 = eqn.sweep(c_var, dt=dt, solver=solver)
    #record the volume integral of the free energy 
    # equivalent to the average value of the free energy for any cell,
    # multiplied by the number of cells and the area of each cell
    # (since this is a 2D domain)
       
    
    for sweeps in range(total_sweeps):
        res = eqn.sweep(c_var, dt=dt, solver=solver)
                    
    if res < res0 * tolerance:
        #print steps
        #print elapsed   
        # anything in this loop will only be executed every 10 steps
        if (steps%10==0):
            print "Saving data"
            if (problem == 'a') or (problem == 'b') or (problem == 'c'):
                save_data(elapsed, c_var, f(c_var).cellVolumeAverage*mesh.numberOfCells*dx*dx, steps)
            
            elif (problem == 'd'):
                #print "saving for d!!"
                save_data(elapsed, c_var, f(c_var).cellVolumeAverage * mesh.cellVolumes.sum()/(.01*r), steps)

        #print "f", f(c_var).cellVolumeAverage * mesh.cellVolumes.sum()
        steps += 1
        elapsed += dt
        dt *= 1.1
        dt = min(dt, dt_max)
        c_var.updateOld()
    else:
        dt *= 0.8
        c_var[:] = c_var.old

# simulation ends
print 'steps reached:', steps


