
# coding: utf-8

# In[ ]:
import numpy as np
import sympy
import fipy as fp
import matplotlib.pyplot as plt
import os
import sys
import json
from fipy.solvers.pysparse import LinearLUSolver as Solver
import resource

memory1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss #get initial memory peak

marker = .5 #add this as a parameter at some point

#load in the json parameter file here
jsonfile = sys.argv[1]

if jsonfile:
    with open(jsonfile, 'rb') as ff:
        params = json.load(ff)
        
else:
    params = dict()
    
print 'my params:', params

#extract the parameters
N = params.get('N', 20)
dx = params.get('dx', 1)  
total_steps = params.get('steps', 2)
sumatra_label = params.get('sumatra_label', '')
total_sweeps = params.get ('sweeps', 2)

c, rho_s, c_alpha, c_beta = sympy.symbols("c_var rho_s c_alpha c_beta")
f_0 = rho_s * (c - c_alpha)**2 * (c_beta - c)**2

sympy.diff(f_0, c, 2)

# command format:
# python 1a_commandline.py 100 2.0
# where nx = 100
# and dx = 2.0
# nx = int(sys.argv[1])
# dx = float(sys.argv[2])

mesh = fp.PeriodicGrid2D(nx=N, ny=N, dx=dx, dy=dx)

c_alpha = 0.3
c_beta = 0.7
kappa = 2.0
M = 5.0
c_0 = 0.5
epsilon = 0.01
rho_s = 5.0
filepath = os.path.join('/data/aem1/new1a/surf-research/Andrey/Data', sumatra_label) # 

# solution variable
c_var = fp.CellVariable(mesh=mesh, name=r"$c$", hasOld=True)

# array of sample c-values: used in f versus c plot
vals = np.linspace(-.1, 1.1, 1000)

x , y = np.array(mesh.x), np.array(mesh.y)

# initial value for square and T domains
c_var[:] = c_0 + epsilon * (np.cos(0.105 * x) * np.cos(0.11 * y) + (np.cos(0.13 * x) * np.cos(0.087 * y))**2 + np.cos(0.025 * x - 0.15 * y) * np.cos(0.07 * x - 0.02 * y))

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
f_data = []
time_data = []
cvar_data = []
step_data = []


def save_data(filename, f, time, cvar, steps, N):
    np.savez(os.path.join(filepath, filename), cvar = cvar.value, time = time, f = f.value, steps = steps, N = N)

# solver equation    
eqn = fp.TransientTerm(coeff=1.) == fp.DiffusionTerm(M * f_0_var(c_var)) - fp.DiffusionTerm((M, kappa))

elapsed = 0.0
steps = 0
dt = 0.01 
tolerance = 1e-1

# controls on how long the simulation runs: steps, duration, or both
duration = 300

c_var.updateOld()
from fipy.solvers.pysparse import LinearLUSolver as Solver
from numpy.testing import assert_almost_equal as npt 
solver = Solver()
print "Starting Solver."
while steps <= total_steps:
    res0 = eqn.sweep(c_var, dt=dt, solver=solver)

    for sweeps in range(total_sweeps):
        res = eqn.sweep(c_var, dt=dt, solver=solver)

        if (elapsed + dt) > marker:
            print 'triggered dt = ',dt
            dt_0 = dt
            dt = marker - elapsed
        if elapsed == marker: ##NEED A BETTER CHECK, CANT CHECK FLOATS
            markerfile = "1a{0}step{1}.mpz".format(N,steps)
            save_data(markerfile, f(c_var).cellVolumeAverage*mesh.numberOfCells*(dx**2), elapsed, c_var, steps, N)
            dt = dt_0
            print 'marker time = ', elapsed
            
    if res < res0 * tolerance:
        
        if ((steps%10)==0):
            print steps
            print elapsed
            
            # record the volume integral of the free energy 
            # equivalent to the average value of the free energy for any cell,
            # multiplied by the number of cells and the area of each cell
            # (since this is a 2D domain)
            file_name = "1a{0}x{1}step{2}".format(N, dx, steps)
            save_data(file_name, f(c_var).cellVolumeAverage*mesh.numberOfCells*(dx**2), elapsed, c_var, steps, N)
             
        steps += 1
        elapsed += dt
        dt *= 1.1
        c_var.updateOld()
    else:
        dt *= 0.8
        c_var[:] = c_var.old

#Keep track of how much memory was used and dump into a txt file
memory2 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss #final memory peak
memory_diff = (memory2 - memory1,) #difference between memory peaks
filename2 = 'memory_usage.txt' #will always be the last file save
np.savetxt(os.path.join(filepath, filename2), memory_diff )


#create 400 mesh to interpolate onto
# m = fp.Grid2D(nx=400, ny=400, dx=.5, dy=.5)
# v = fp.CellVariable(mesh=m)

# #interpolate the latest c_var values (hopefully at equillibrium) to the empty 400 grid
# v = c_var([mesh.x, mesh.y], order-1)

#save the interpolated values to a file
# interp_name = "1a{0}x{1}step{2}Interp400".format(N, dx, steps)
# np.savez(os.path.join(filepath, interp_name), cvar = v)

# simulation ends
print 'elapsed_time:', elapsed

