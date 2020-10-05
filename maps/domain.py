"""Domain module

This module contains classes to represent the integration domain, the pops and the resources
"""

import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.ndimage import gaussian_filter
from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image


def restrict(In, dx, new_dx):
    """To project on a coarser grid"""
    step = int(new_dx/dx)
    return In[::step, ::step]


class Domain:
    """The domain of integration"""

    def __init__(self, name):
        """Constructor of the domain

        name -- a string with the name of the domain files
        """

        f = open(name+'/specs.json',)
        specs = json.load(f)

        # Reading all the images
        self.I_start = plt.imread(name+'/bound.png')
        self.I_topo_start = (1-plt.imread(name+'/topo.png'))
        terrain = (plt.imread(name+'/terrain.png')*255).astype(np.uint8)

        # Building the ressource base regeneration map
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
        """To project on a coarser grid

        new_dx -- The new space-step (km)
        """

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
    """Population class"""

    def __init__(self, start_loc, start_number, color, gamma, D0, KN, c, drift0, shape, area, dx):
        """Population constructor

        start_loc -- array with starting location coordonates
        start_number -- initial population density (pop/km^2)
        color -- array with the color of the pop on the map (RGB)
        gamma -- death rate density (%/years)
        D0 -- diffusion coefficient in the best conditions (km^2/years)
        KN -- carrying capacity in the best conditions (pop/km)
        c -- resource transformation capacity ( pop / res)
        drift0 -- base resource attraction (km^2 / res*year)
        shape -- dimensions of the domain (int, int)
        area -- area of the domain (km^2)
        dx -- space step of the domain (km)
        """

        self.pi = np.zeros(shape)

        for i in np.arange(-10,10):
            for j in np.arange(-10,10):
                self.pi[start_loc[0]+i, start_loc[1]+j] = start_number / area

        self.gamma = gamma*np.ones(shape)
        self.D = D0*np.ones(shape)
        self.KN = KN*np.ones(shape) / area
        self.repro = np.zeros(shape)
        self.color = color
        self.dx = dx
        self.v = 0.002
        self.c = c
        self.drift0 = drift0

    def mask(self):
        """returns a boolean mask of the occupation zone"""

        out = np.zeros_like(self.pi)
        out[np.where(self.pi == 0)] = 1
        return True*out

    def area(self):
        """returns the area occupied in km^2"""

        return np.sum(self.mask())*(self.dx**2)

    def tot(self):
        """returns the total population in pop"""
        return np.sum(self.pi)*self.area()

    def resize(self, new_dx):
        """projection on a coarser grid"""

        if new_dx > self.dx:
            self.pi = restrict(self.pi, self.dx, new_dx)
            self.gamma = restrict(self.gamma, self.dx, new_dx)
            self.D = restrict(self.D, self.dx, new_dx)
            self.KN = restrict(self.KN, self.dx, new_dx)
            self.repro = restrict(self.repro, self.dx, new_dx)
            self.dx = new_dx
        else:
            print("only coarser grid are authorized !")

    def alpha(self, u):
        return np.exp(-(self.v-u)**2/4)

    def G(self, a, U, RHO, PI, i):
        conso = self.c*np.sum([a[i,j]*RHO[j].rho for j in range(RHO.shape[0])], axis=0)
        death = -self.gamma
        war = -np.sum([self.alpha(U[i])*PI[i].pi for i in range(PI.shape[0])], axis=0)
        return conso+death+war

    def DN(self, r):
        """Diffusion coefficient (see Gorban and all)"""
        return self.D*np.exp(-(r/np.linalg.norm(r)))

    def colorize(self):
        """Returns the colormap of occupied zones"""

        out = self.color*np.ones((self.pi.shape[0], self.pi.shape[1], 3))
        #fact = (1-np.exp(-self.pi/np.max(self.pi)))
        fact = np.zeros_like(self.pi)
        fact[np.where(self.pi>0.1*np.max(self.pi))]=1

        out[:,:,0] *= fact
        out[:,:,1] *= fact
        out[:,:,2] *= fact
        return out


class Res:
    """Resources class"""

    def __init__(self, r0, KR, shape, area):
        """Constructor of the class

        r0 -- The renewal rate of resources in best conditions (%/year)
        KR -- The carrying capacity density in best conditions (res/km)
        shape -- dimensions of the domain (int, int)
        area -- area of the domain (km^2)
        """

        self.rho = (KR/area)*np.ones(shape)
        self.r = r0*np.ones(shape)
        self.KR = KR*np.ones(shape)

    def tot(self, area):
        """returns the total amount of resource"""

        return np.sum(self.rho)*area

    def resize(self, new_dx, old_dx):
        """projection on a coarser grid"""

        if new_dx > old_dx:
            self.rho = restrict(self.rho, old_dx, new_dx)
            self.r = restrict(self.r, old_dx, new_dx)
            self.KR = restrict(self.KR, old_dx, new_dx)
        else:
            print("only coarser grid are authorized !")
