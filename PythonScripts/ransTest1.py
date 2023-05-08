#from read_unstructured_grid_2 import get_geo
from lib2to3.pgen2.token import LESSEQUAL
from vtklb import vtklb
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
#import vtk_tools

def interface_profile(centerpos, pos, WInv):
    return 0.5*(1-np.tanh(2*WInv*(pos-centerpos)))

def cirular_shift(geo, dr):
    M = np.roll(geo, dr[0], axis=0)
    M = np.roll(M, dr[1], axis=1)
    return M  

def figPlot2D(figno,fieldX):
    plt.figure(figno)
    ax = plt.pcolormesh(fieldX.transpose())
    plt.axis('equal')
    plt.axis('off')
    plt.colorbar()

plt.rcParams['text.usetex'] = True

np.random.seed(0)
# system size: nx  ny
system_size = (300, 50)
# set the size of the geometri
geo = np.ones(system_size, dtype=int)
# partition the system in two
nproc = 4
for x in range(nproc):
    print(str(x))
    y = x + 1
    geo[int((y-1)*system_size[0]/nproc):int(y*system_size[0]/nproc), :] = y

Cmu= 0.09
l_turb = 10

kInit = 1e-4 #1e-4 #0.2
EInit = Cmu**0.75 * kInit**1.5 / l_turb

print('kInit = '+str(kInit))
print('EInit = '+str(EInit))

# path to your badchimp folder
path_badchimp = "/home/AD.NORCERESEARCH.NO/olau/Programs/Git/BADChIMP-cpp/"



#geo[:, 0] = 0
#geo[:, -1] = 0
x0 = int(0.25*system_size[0])

y0 = int(0.5*system_size[1])
RGeo = 6 #int(0.5*system_size[0])


for x in range(system_size[0]):
    for y in range(system_size[1]):
        rSq = (x-x0)*(x-x0) + (y-y0)*(y-y0)
        r = np.sqrt(rSq)
        if r <= RGeo :
            geo[x,y] = 0



#geo[:15,:20] = 0

vtk = vtklb(geo, "D2Q9", "xy", "tmp", path_badchimp + "input/mpi/") 

rho = np.ones(system_size)
rhoK = kInit*np.ones(system_size)
rhoE = EInit*np.ones(system_size)




vtk.append_data_set("init_rho", rho)
vtk.append_data_set("init_rhoK", rhoK)
vtk.append_data_set("init_rhoEpsilon", rhoE)


figPlot2D(3, geo)

rho[geo==0] = np.nan
figPlot2D(4, rho)

rhoK[geo==0] = np.nan
figPlot2D(5, rhoK)

rhoE[geo==0] = np.nan
figPlot2D(6, rhoE)




plt.show()