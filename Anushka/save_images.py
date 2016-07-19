import matplotlib.pyplot as plt
import numpy as np
import fipy as fp
import glob

nx = 200
dx = 1.0

#mesh = fp.PeriodicGrid2D(nx=nx, ny=nx, dx=dx, dy=dx)
#mesh = fp.Grid2D(nx=nx, ny=nx, dx=dx, dy=dx)
mesh = fp.Grid2D(Lx=20., Ly=100.0, nx=nx / 5, ny=nx) + (fp.Grid2D(Ly=20.0, Lx=100.0, nx=nx, ny=nx / 5) + [[-40],[100]])

# load and cut off values at given duration
file_dir = '/data/and9/surf-research/Anushka/1c'
filename = file_dir + '/1c200*.npz'
newest = max(glob.iglob(filename), key=os.path.getctime)


times = np.load(newest)['time']
c_np = np.load(newest)['c_var']
f = np.load(newest)['f']

c_np_0 = c_np[0]

c_var = fp.CellVariable(value = c_np_0, mesh = mesh)
#viewer = fp.Viewer(c_var, title = "Time")


fig = plt.figure(figsize = (10, 5))
#fig.suptitle('Periodic Boundary Conditions on a Square Domain', fontsize = 15)
#fig.suptitle('Fixed-Flux Boundary Conditions on a Square Domain', fontsize = 15)
fig.suptitle('Fixed-Flux Boundary Conditions on a T-Shaped Domain', fontsize = 15)

axes1 = fig.add_subplot(121)
viewer = fp.Viewer(c_var, axes=axes1, colorbar=None, title = 'Concentration Distribution', fontsize = 10)
viewer.axes.get_xaxis().set_visible(False)
viewer.axes.get_yaxis().set_visible(False)

axes2 = fig.add_subplot(122)
axes2.set_xlabel('Time', fontsize = 10)
axes2.set_ylabel('Free Energy', fontsize = 10)
axes2.set_title('Total Free Energy Evolution', fontsize = 10)
axes2.set_yscale('symlog')
axes2.set_xscale('symlog')
plt.subplots_adjust(left=.05, right=.95, top=.85, bottom=.1, wspace = .2)
    
for index in range(len(times)):
    
    current_time = times[index]
#    print current_time
    
    f_current = f[:(index+1)]
    times_current = times[:(index+1)]
    c_np_current = c_np[index]
    c_var[:] = c_np_current
#    c_var = fp.CellVariable(value = c_np_current, mesh = mesh)

    
    viewer.plot()

    axes2 = fig.add_subplot(122)
    axes2.set_xlim(0, times[-1] + 10)
    axes2.set_ylim(0, f[0] + 10)
 #   print times_current
#    print f_current

    axes2.plot(times_current, f_current, color = 'b')


#    fig.show() 
 
#    plt.savefig('1b_images/image{0}.png'.format(str(index).rjust(5, '0')))
    plt.savefig('/tmp/1c_images/image{0}.png'.format(str(index).rjust(5, '0')))
    print "saving ", 'image{0}.png'.format(str(index).rjust(5, '0')) 

