"""Grid module

This module combines domain, pops and resources to make the simulation
"""


import numpy as np


class Grid:
    """The grid class"""

    def __init__(self,  PI, RHO, Domain, dx=None):
        """Constructor of the grid

        PI -- an array of pops
        RHO -- an array of resources
        Domain -- The domain on which the simulation will be made
        dx -- The sharpness of the discretization
        """

        self.RHO = RHO
        self.PI = PI
        self.dom = Domain

        for pi in self.PI:
            pi.D *= self.dom.I_topo
            pi.KN *= self.dom.I_r

        for rho in self.RHO:
            rho.KR *= self.dom.I_topo
            rho.r *= self.dom.I_r

        if dx != None:
            self.dom.resize(dx)
            for pi in self.PI:
                pi.resize(dx, self.dom.dx_start)
            for rho in self.RHO:
                rho.resize(dx, self.dom.dx_start)

        self.a = np.zeros((self.PI.shape[0], self.RHO.shape[0], self.dom.shape[0], self.dom.shape[1]))

        self.time = 0

    def update(self):
        """Update the grid by one time step"""

        self.time += 1
        fluctuations = 2*(1-np.sin(100*self.time))

        for i,pi in enumerate(self.PI):
            pi.update(self.a, self.RHO, self.PI,i)
            pi.pi *= self.dom.I

        for j,rho in enumerate(self.RHO):
            rho.update(fluctuations, self.a, self.RHO, self.PI,j)
            rho.rho *= self.dom.I

        for i in range(self.PI.shape[0]):
            for j in range(self.RHO.shape[0]):
                self.a[i,j] += self.RHO[j].rho
                self.a[i,j] /= np.max(self.a[i,j])
                self.a[i,j] *= 0.0002*self.dom.I

    def get_img(self):
        """returns the repartition of all pops"""

        out = np.ones((self.dom.shape[0], self.dom.shape[1],3))
        out[:,:,0] = self.dom.I
        out[:,:,1] = self.dom.I
        out[:,:,2] = self.dom.I

        for pi in self.PI:
            out += pi.colorize()
        return (out*255).astype(np.uint8)
