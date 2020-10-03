import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.ndimage import gaussian_filter
from maps.mathutils import grad_grad, lap


from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image

def rgb_to_bw(Im):
    return 0.2989 * Im[:,:,0] + 0.5870 * Im[:,:,1] + 0.1140 * Im[:,:,2]

def restrict(In, dx, new_dx):
    step = int(new_dx/dx)
    return In[::step, ::step]


class Domain:
    def __init__(self, name):

        f = open(name+'/specs.json',)
        specs = json.load(f)

        self.I_start = plt.imread(name+'/bound.png')
        self.I_topo_start = (1-plt.imread(name+'/topo.png'))
        terrain = (plt.imread(name+'/terrain.png')*255).astype(np.uint8)

        self.I_r_start = np.zeros(self.I_start.shape)

        for t in specs["terrain"]:
            self.I_r_start[np.where(np.all(terrain == np.array(t["color"]), axis=-1))] = t["prosperity"]


        self.dx_start = specs["dx"]

        self.I = self.I_start
        self.I_topo = self.I_topo_start
        self.I_r = self.I_r_start

        self.dx = specs["dx"] #km
        self.area = np.sum(self.I)*(self.dx**2)
        self.shape = self.I.shape

    def resize(self, new_dx):
        if new_dx > self.dx_start:
            self.I = restrict(self.I_start, self.dx_start, new_dx)
            self.I_topo = restrict(self.I_topo_start, self.dx_start, new_dx)
            self.I_r = restrict(self.I_r_start, self.dx_start, new_dx)

            self.dx = new_dx
            self.area = np.sum(self.I)*(self.dx**2)
            self.shape = self.I.shape
        else:
            print("only coarser grid are authorized !")



class Pop:
    def __init__(self, start_loc, start_number, color, gamma, D0, KN, shape, area):
        self.pi = np.zeros(shape)

        for i in np.arange(-10,10):
            for j in np.arange(-10,10):
                self.pi[start_loc[0]+i, start_loc[1]+j] = start_number / area

        self.gamma = gamma*np.ones(shape)
        self.D = D0*np.ones(shape)
        self.KN = KN*np.ones(shape) / area
        self.repro = np.zeros(shape)
        self.color = color

    def mask(self):
        out = np.zeros_like(self.pi)
        out[np.where(self.pi == 0)] = 1
        return True*out

    def area(self, dx):
        return np.sum(self.mask())*(dx**2)

    def tot(self, area):
        return np.sum(self.pi)*area

    def resize(self, new_dx, old_dx):
        if new_dx > old_dx:
            self.pi = restrict(self.pi, old_dx, new_dx)
            self.gamma = restrict(self.gamma, old_dx, new_dx)
            self.D = restrict(self.D, old_dx, new_dx)
            self.KN = restrict(self.KN, old_dx, new_dx)
            self.repro = restrict(self.repro, old_dx, new_dx)
        else:
            print("only coarser grid are authorized !")

    def update(self,a,RHO,PI,i):

        def DN(r):
            return self.D*np.exp(-(r/np.linalg.norm(r)))

        conso = np.sum([a[i,j]*RHO[j].rho for j in range(RHO.shape[0])], axis=0)
        death = -self.gamma
        self.repro = conso+death

        migration =  DN(self.repro)*lap(self.pi)
        shift = 0.000001*self.pi*lap(np.sum([rho.rho for rho in RHO], axis=0))

        self.pi += self.pi*self.repro + migration + shift

    def colorize(self):
        out = self.color*np.ones((self.pi.shape[0], self.pi.shape[1], 3))
        #fact = (1-np.exp(-self.pi/np.max(self.pi)))
        fact = np.zeros_like(self.pi)
        fact[np.where(self.pi>0.1*np.max(self.pi))]=1

        out[:,:,0] *= fact
        out[:,:,1] *= fact
        out[:,:,2] *= fact
        return out


class Res:
    def __init__(self, r0, KR, shape, area):
        self.rho = (KR/area)*np.ones(shape)
        self.r = r0*np.ones(shape)
        self.KR = KR*np.ones(shape)

    def tot(self, area):
        return np.sum(self.rho)*area

    def resize(self, new_dx, old_dx):
        if new_dx > old_dx:
            self.rho = restrict(self.rho, old_dx, new_dx)
            self.r = restrict(self.r, old_dx, new_dx)
            self.KR = restrict(self.KR, old_dx, new_dx)
        else:
            print("only coarser grid are authorized !")

    def update(self,fluctuations,a,RHO,PI,j):
        renew = self.r
        thresh = -(self.r*self.rho)/self.KR
        conso = -np.sum([a[k,j]*PI[k].pi for k in range(PI.shape[0])], axis=0)

        repro_res = (fluctuations*renew+thresh+conso)

        self.rho += self.rho*repro_res
