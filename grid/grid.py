"""Grid module

This module combines domain, pops and resources to make the simulation
"""


import numpy as np
from grid.mathutils import grad_grad, lap


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
                pi.resize(dx)
            for rho in self.RHO:
                rho.resize(dx, self.dom.dx_start)

        self.a = np.zeros((self.PI.shape[0], self.RHO.shape[0], self.dom.shape[0], self.dom.shape[1]))

        self.time = 0

    def update(self):
        """Update the grid by one time step"""

        self.time += 1
        fluctuations = 2*(1-np.sin(100*self.time))

        U = [pi.v for pi in self.PI]

        for i,pi in enumerate(self.PI):

            pi.repro = pi.G(self.a, U,self.RHO,self.PI, i)
            migration = pi.DN(pi.repro)*lap(pi.pi, pi.dx)
            shift = 0.000001*pi.pi*lap(np.sum([rho.rho for rho in self.RHO], axis=0), pi.dx)

            pi.pi += pi.pi*pi.repro + migration + shift
            pi.pi *= self.dom.I

            pi.v += np.sum([np.sqrt(2)*((pi.v-U[i])/2)*pi.alpha(U[i]) for i in range(self.PI.shape[0])])

        for j,rho in enumerate(self.RHO):
            renew = rho.r
            thresh = -(rho.r*rho.rho)/rho.KR
            conso = -np.sum([self.a[k,j]*self.PI[k].pi for k in range(self.PI.shape[0])], axis=0)

            repro_res = (fluctuations*renew+thresh+conso)

            rho.rho += rho.rho*repro_res
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
