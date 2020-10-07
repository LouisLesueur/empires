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

        self.time = 0

    def update(self):
        """Update the grid by one time step"""

        self.time += 1
        fluctuations = 2*(1-np.sin(100*self.time))

        U = [pi.v for pi in self.PI]

        for i,pi in enumerate(self.PI):

            pi.repro = pi.G(U,self.RHO,self.PI)
            migration = pi.DN(pi.repro)*lap(pi.pi, pi.dx)
            shift = pi.drift0*pi.pi*lap(np.sum([rho.rho for rho in self.RHO], axis=0), pi.dx)

            pi.pi += pi.pi*pi.repro + migration + shift

            pi.v += pi.partG(U,self.RHO,self.PI)

        for j,rho in enumerate(self.RHO):
            renew = rho.r
            thresh = -(rho.r*rho.rho)/rho.KR
            conso = -np.sum([rho.prop*self.PI[k].pi for k in range(self.PI.shape[0])], axis=0)

            repro_res = (fluctuations*renew+thresh+conso)

            rho.rho += rho.rho*repro_res


    def get_img(self):
        """returns the repartition of all pops"""

        out = np.ones((self.dom.shape[0], self.dom.shape[1],3))

        out[:,:,0] = self.dom.I
        out[:,:,1] = self.dom.I
        out[:,:,2] = self.dom.I

        for pi in self.PI:
            out += pi.colorize()

        return (out*255).astype(np.uint8)
