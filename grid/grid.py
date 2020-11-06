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
        self.D = np.zeros_like(self.dom.I)
        self.P = np.zeros_like(self.dom.I)
        self.Z = np.zeros_like(self.dom.I)
        self.G = np.zeros_like(self.dom.I)
        self.S = np.zeros_like(self.dom.I)
        self.C = np.zeros((len(Ns.start_loc),*self.dom.I.shape))
        self.R = 0.75*Rs.Rmax*np.ones_like(self.dom.I)*self.dom.I
        self.Idx = np.zeros_like(self.dom.I)-1



        self.x = np.zeros((*self.dom.I.shape,2))
        for i in range(self.dom.I.shape[0]):
            for j in range(self.dom.shape[1]):
                self.x[i,j] = np.array([i,j])*self.dx

        self.states = []
        self.time = 0
        self.dt = dt

        self.r0 = Rs.r0*self.dom.I_r
        self.Rmax = 2*Rs.Rmax*self.dom.I_topo

        self.c = 1
        self.eps = 1
        self.c1 = 3
        self.alpha = 0.2
        self.sig = 100
        self.h = 100
        self.z = 100

        for i,idx in enumerate(Ns.start_loc):
            if self.dom.I[idx[0],idx[1]]:
                self.C[i,idx[0],idx[1]] = 1
                self.N[idx[0],idx[1]] = Ns.N_start[i]
                self.Idx[idx[0],idx[1]] = i
                for id in range(-1,2):
                    for jd in range(-1,2):
                        if self.dom.I[idx[0]+id, idx[1]+jd]:
                            self.Idx[idx[0]+id, idx[1]+jd] = i
                            self.N[idx[0]+id, idx[1]+jd] = Ns.N_start[i]

            self.states.append(State(np.random.rand(3),i,self.c,
                                     self.eps, self.c1,self.alpha,
                                     self.sig, self.h,self.z))
        self.N*=self.dom.I


        self.c0 = Ns.c0
        self.Rdem = Ns.Rdem
        self.Zdem = self.c*Ns.Rdem
        self.n0 = Ns.n0
        self.k0 = Ns.k0
        self.Nbar = Ns.Nbar
        self.gamma = Ns.gamma
        self.bary_map = np.zeros((*self.dom.I.shape,2))
        self.areas = np.array([np.sum((self.Idx == s.idx).astype(np.int32))*(self.dx**2) for s in self.states])
        self.barycenters = np.array([np.sum(self.x[np.where(self.Idx==s.idx)], axis=0)/len(self.x[np.where(self.Idx==s.idx)]) for s in self.states])

        for s in self.states:
            self.bary_map[np.where(self.Idx == s.idx)] = self.barycenters[s.idx]


    def update(self):
        """Update the grid by one time step"""
        self.time += self.dt
        I_filter = self.Idx>-1

        # Renewal
        self.R += self.r0*self.R*(1-(self.R/(self.Rmax+1e-6)))*self.dt

        # Production
        self.Z += self.c*self.eps*self.R*self.dt*I_filter
        self.R -= self.eps*self.R*self.dt*I_filter

        #Taxes
        self.Zpub = self.alpha*self.Z*self.dt*I_filter
        self.Zpriv = (1-self.alpha)*self.Z*self.dt*I_filter

        #Consumption
        self.Zpub -= ((self.c0*self.R)/(self.Rdem + self.R))*self.dt*I_filter
        satisfaction = (self.Zpub - self.Zdem)*I_filter

        #Cultural assimilation
        self.C += self.c1* np.einsum('ijk,ijk->ijk',self.C,(1-self.C))*I_filter

        #Demography
        k = (1 + self.Zpub/(self.Zdem + self.Zpub))*self.Nbar*I_filter
        self.G = self.n0*(1-(self.N/(k+1e-6)))
        self.N += self.G*self.N*self.dt*I_filter

        #Migration
        self.D = 0.1
        self.N += self.D*lap(self.N, self.dx)*self.dt
        for s in self.states:
            z = np.zeros_like(self.Idx, dtype=np.bool)
            z[np.where(self.Idx == s.idx)] = True
            z = (self.D*lap(z, self.dx)).astype(np.bool)
            z[np.where(self.N>0)] = False
            self.Idx[np.where(z==True)] = s.idx


        #Geographical update
        self.areas = np.array([np.sum((self.Idx == s.idx).astype(np.int32))*(self.dx**2) for s in self.states])
        self.barycenters = np.array([np.sum(self.x[np.where(self.Idx==s.idx)], axis=0)/(len(self.x[np.where(self.Idx==s.idx)]+1e-6)) for s in self.states])
        for s in self.states:
            self.bary_map[np.where(self.Idx == s.idx)] = self.barycenters[s.idx]

        #Asalyia
        # self.S += self.S*(1-self.S)*np.exp(np.linalg.norm(self.x-self.bary_map)**2 / self.sig**2)*np.exp((satisfaction)**2 / (1.5*self.Zdem)**2)
        # for id in np.where(self.S < 0.2):
        #     self.Idx[id] = len(self.states)
        #     self.states.append(State(np.random.rand(3),len(self.states),self.c,
        #                              self.eps, self.c1,self.alpha,
        #                              self.sig, self.h,self.z))

        #Power
        self.P = np.exp(-np.linalg.norm(self.x-self.bary_map)/self.h)*(self.Zpub / self.Zdem)*(self.dx**2)*self.dt


    def get_img(self):
        """returns the repartition of all pops"""

        out = np.ones((*self.dom.I.shape,3))

        out[:,:,0] = self.dom.I
        out[:,:,1] = self.dom.I
        out[:,:,2] = self.dom.I

        for s in self.states:
            out[np.where(self.Idx==s.idx)] = s.color
        #out[np.where(self.N>0)] = (self.N/np.max(self.N))*np.array([1,0,0])


        return (out*255).astype(np.uint8)
        #return self.D
