import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import fipy as fp
from fipy import numerix as nmx

r = 100
numCellsDesired = 2500

epsilon = 0.05
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


# load and cut off values at given duration
times = np.load('/data/and9/surf-research/Anushka/1d/1d50.npz')['time']
c_np = np.load('/data/and9/surf-research/Anushka/1d/1d50.npz')['c_var']
f = np.load('/data/and9/surf-research/Anushka/1d/1d50.npz')['f']

c_np_0 = c_np[0]

c_var = fp.CellVariable(value = c_np_0, mesh = mesh)
viewer = fp.Viewer(c_var, title = "Time")

fig = plt.figure(figsize = (10, 5))
fig.suptitle('Spherical Domain', fontsize = 15)

axes1 = fig.add_subplot(121)

viewer = fp.Viewer(c_var, colorbar=None, title = 'Concentration Distribution', fontsize = 10) 
axes1.axes.get_xaxis().set_visible(False)
axes1.axes.get_yaxis().set_visible(False)

axes2 = fig.add_subplot(122)
axes2.set_xlabel('Time', fontsize = 10)
axes2.set_ylabel('Free Energy', fontsize = 10)
axes2.set_title('Total Free Energy Evolution', fontsize = 10)
axes2.set_yscale('symlog')
axes2.set_xscale('symlog')
plt.subplots_adjust(left=.05, right=.95, top=.85, bottom=.1, wspace = .2)

for index in range(1):
    
    current_time = times[index]  
    f_current = f[:(index+1)]
    times_current = times[:(index+1)]
    c_np_current = c_np[index]
    c_var[:] = c_np_current

    viewer.plot('viewer.png')
    img = mpimg.imread('viewer.png')
    axes1.imshow(img[: , 160:385])
    
    axes2 = fig.add_subplot(122)
    axes2.set_xlim(0, times[-1] + 10)
    axes2.set_ylim(0, f[0] + 10)
    axes2.plot(times_current, f_current, color = 'b')
        
    fig.savefig('/tmp/1d_images/image{0}.png'.format(str(index).rjust(5, '0')))
   
