"""Grid module

This module combines domain, pops and resources to make the simulation
"""


import numpy as np
from grid.mathutils import grad_grad, lap
from maps.domain import State

class Grid:
    """The grid class"""

    def __init__(self,  Ns, Rs, Domain, dt):
        """Constructor of the grid

        N -- pops
        Rs -- resources
        Domain -- The domain on which the simulation will be made
        dt -- Time step (year)
        """

        self.dom = Domain
        self.dx = Domain.dx

        self.N = np.zeros_like(self.dom.I)
        self.P = np.zeros_like(self.dom.I)
        self.Z = np.zeros_like(self.dom.I)
        self.C = np.zeros((len(Ns.start_loc),*self.dom.I.shape))
        self.R = 0.75*Rs.Rmax*np.zeros_like(self.dom.I)
        self.Idx = np.zeros_like(self.dom.I)-1

        self.x = np.array([[[i,j] for i in range(self.dom.I.shape[0])] for j in range(self.dom.I.shape[1])])*self.dx

        self.states = []
        self.time = 0
        self.dt = dt

        self.r0 = Rs.r0*self.dom.I_topo
        self.Rmax = Rs.Rmax*self.dom.I_r

        self.c = 1
        self.eps = 0.7
        self.c1 = 3
        self.alpha = 0.2
        self.sig = 100
        self.h = 100
        self.z = 100

        for i,idx in enumerate(Ns.start_loc):
            self.C[i,idx] = 1
            self.N[idx] = Ns.N_start[i]
            self.Idx[idx] = i
            self.states.append(State(np.random.rand(3),i,self.c,
                                     self.eps, self.c1,self.alpha,
                                     self.sig, self.h,self.z))

        self.c0 = Ns.c0
        self.Zdem = self.c*Ns.Rdem
        self.n0 = Ns.n0
        self.k0 = Ns.k0
        self.Nbar = Ns.Nbar
        self.gamma = Ns.gamma
        self.barycenters = np.zeros((*self.dom.I.shape,2))
        self.areas = np.zeros(len(self.states))


    def update(self):
        """Update the grid by one time step"""
        self.time += self.dt

        # Renewal
        self.R += self.r0*self.R*(1-(self.R/self.Rmax))*self.dt

        # Production
        self.Z += self.c*self.eps*self.R*self.dt
        self.R -= self.eps*self.R*self.dt

        #Taxes
        self.Zpub = self.alpha*self.Z*self.dt
        self.Zpriv = (1-self.alpha)*self.Z*self.dt

        #Consumption
        self.Zpub -= ((self.c0*self.R)/(self.Rdem + self.R))*self.dt
        satisfaction = self.Zpub - self.Zdem

        #Cultural assimilation
        self.C += self.c1* np.einsum('ijk,ijk->ijk',self.C,(1-self.C))

        #Demography
        k = self.k0*(1 + self.Zpub/(self.Zdem + self.Zpub))*self.Nbar
        G = self.n0*(1-(self.N/k))
        self.N += G*self.N*self.dt

        #Migration
        D = 0.5*np.tanh(self.gamma*self.N)
        self.N += D*lap(self.N)*(self.Idx == -1)

        #Asalyia
        bary = np.array([np.sum(np.where(A==s.idx)*self.dx,axis=0)/len(np.where(A==s.idx))
                                     for s in self.states])
        for s in self.states:
            self.barycenters[np.where(self.Idx == s.idx),:] = bary[s.idx]


        self.S += self.S*(1-self.S)*np.exp((self.x-self.barycenters)**2 / self.sig**2)*np.exp((satisfaction)**2 / (1.5*self.Zdem)**2)
        for id in np.where(self.S < 0.2):
            self.Idx[id] = len(self.states)
            self.states.append(State(np.random.rand(3),len(self.states),self.c,
                                     self.eps, self.c1,self.alpha,
                                     self.sig, self.h,self.z))

        #Power
        self.areas = np.array([np.sum(np.where(A==s.idx)*self.dx**2)
                                     for s in self.states])
        self.P = np.einsum('i,ijk -> ijk',self.areas,self.S)*np.exp(-np.linalg.norm(self.x-self.barycenters)/self.h)*(self.Zpub / self.Zdem)




    def get_img(self):
        """returns the repartition of all pops"""

        out = np.ones((*self.dom.I.shape,3))

        out[:,:,0] = self.dom.I
        out[:,:,1] = self.dom.I
        out[:,:,2] = self.dom.I

        for s in self.states:
            out[np.where(self.Idx==s.idx),0] = s.color[0]
            out[np.where(self.Idx==s.idx),1] = s.color[1]
            out[np.where(self.Idx==s.idx),2] = s.color[2]

        return (out*255).astype(np.uint8)
