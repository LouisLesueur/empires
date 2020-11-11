"""Grid module

This module combines domain, pops and resources to make the simulation
"""


import numpy as np
from grid.mathutils import grad_grad, lap, compute_bound, bound_lap
from maps.domain import State
from scipy.ndimage import convolve

class Grid:
    """The grid class"""

    def __init__(self,  Ns, Rs, Domain, alpha, Nbar, dt):
        """Constructor of the grid

        N -- pops
        Rs -- resources
        Domain -- The domain on which the simulation will be made
        alpha -- Initial values for each state
        Nbar -- barbarian population (pop)
        dt -- Time step (year)
        """

        self.dom = Domain
        self.dx = Domain.dx

        self.n = self.dom.I.shape[0]
        self.m = self.dom.I.shape[1]

        self.x, self.y = np.meshgrid(np.arange(self.m)*self.dx, np.arange(self.n)*self.dx)

        self.N = np.zeros_like(self.dom.I)
        self.Rpub = np.zeros_like(self.dom.I)
        self.Rpriv = np.zeros_like(self.dom.I)
        self.conso = np.zeros_like(self.dom.I)
        self.satisfaction = np.zeros_like(self.dom.I)
        self.Idx = np.zeros_like(self.dom.I)-1



        self.time = 0

        self.r0 = Rs.r0*self.dom.I_topo*self.dom.I
        self.Rmax = Rs.Rmax*self.dom.I_topo
        self.R = 0.75*self.Rmax*self.dom.I

        self.alpha = alpha
        self.states = {}
        for i,idx in enumerate(Ns.start_loc):
            self.state_idmax = i
            if self.dom.I[idx[0],idx[1]]:
                self.N[idx[0],idx[1]] = Ns.N_start[i]
                self.Idx[idx[0],idx[1]] = i
                for id in range(-1,2):
                    for jd in range(-1,2):
                        if self.dom.I[idx[0]+id, idx[1]+jd]:
                            self.Idx[idx[0]+id, idx[1]+jd] = i
                            self.N[idx[0]+id, idx[1]+jd] = Ns.N_start[i]

            self.states[i] = State(np.random.rand(3),i, self.alpha)

        self.N*=self.dom.I


        self.c0 = Ns.c0
        self.Rdem = Ns.Rdem
        self.n0 = Ns.n0
        self.chi = Ns.chi

        self.D = Ns.D
        self.drift = Ns.drift

        self.dt = dt
        if(self.dt >= (self.dx**2)/(4*self.D)):
            print("Warning, too large dt !")

        self.Nmax = Ns.Nmax*self.dom.I_topo
        self.Nbar = Nbar
        self.bary_map = np.zeros((*self.dom.I.shape,2))
        self.conflicts = np.zeros((*self.dom.I.shape,2))
        for s in self.states.values():
            self.bary_map[np.where(self.Idx == s.idx)] = s.barycenter(self.x, self.Idx)


    def update(self):
        """Update the grid by one time step"""

        self.time += self.dt

        I_filter = self.Idx>-1

        #Geographical update
        #---------------------------------------------------------------------------------------------------------
        remove = []
        for s in self.states.values():
            if s.area(self.Idx, self.dx)==0:
                remove.append(s.idx)
        for r in remove:
            del self.states[r]

        if len(self.states) == 0:
            return

        #Consumption
        #---------------------------------------------------------------------------------------------------------
        self.conso = self.c0*(self.R/(self.Rdem + self.R))*self.N
        self.satisfaction = (self.R - self.Rdem)

        # Renewal
        #---------------------------------------------------------------------------------------------------------
        renew = self.r0
        limit = -self.r0*(self.R/self.Rmax)

        dR = (renew+limit)*self.R-self.conso
        self.R += dR*self.dt


        #Taxes
        #---------------------------------------------------------------------------------------------------------
        self.Rpub = self.alpha*self.R*I_filter
        self.Rpriv = (1-self.alpha)*self.R*I_filter

        self.satisfaction[I_filter] = (self.Rpriv - self.Rdem)[I_filter]

        #Revolts
        #---------------------------------------------------------------------------------------------------------
        insatisfaction = -13
        if len(np.array(np.where(self.satisfaction < insatisfaction)).T) > 0:
            self.Idx[np.where(self.satisfaction < insatisfaction)] = self.state_idmax
            self.states[self.state_idmax]=State(np.random.rand(3),self.state_idmax,self.alpha)
            self.state_idmax += len(np.array(np.where(self.satisfaction < insatisfaction)).T)


        #Demography (Bazykin model + migrations)
        #---------------------------------------------------------------------------------------------------------
        dist = np.array([[np.sqrt(2),1,np.sqrt(2)],
                         [1,         0,         1],
                         [np.sqrt(2),1,np.sqrt(2)]])*(self.dx)

        death = - self.n0
        limit = - (self.N/self.Nmax)
        G = death+limit

        migration = self.D*lap(self.N, self.dx) + self.drift*self.N*lap(self.R, self.dx)
        self.Idx, self.conflicts = bound_lap(self.Idx)
        self.Idx[np.where(self.N<=self.Nbar)] = -1

        dN = (G*self.N +self.chi*self.conso+migration)*self.dt
        self.N += dN*self.dom.I

        #Wars
        #---------------------------------------------------------------------------------------------------------
        common_boundaries = np.array(np.where(self.Idx == -2)).T
        done = []
        if len(common_boundaries) > 0:
            Nwars = 10
            for _ in range(Nwars):
                battle = common_boundaries[np.random.randint(0, len(common_boundaries))]
                if battle not in np.array(done):
                    protagonists = []
                    power = []
                    for i in range(-1,2):
                        for j in range(-1,2):
                            if 0 < battle[0]+i < self.Idx.shape[0] and 0 < battle[1]+j <  self.Idx.shape[1]:
                                if self.Idx[battle[0]+i, battle[1]+j] != -1:
                                    protagonists.append(int(self.Idx[battle[0]+i, battle[1]+j]))
                                    power.append(np.mean(self.Rpub[self.Idx == int(self.Idx[battle[0]+i, battle[1]+j])]))
                    self.Idx[battle[0],battle[1]] = protagonists[np.argmax(power)]
                    done.append(battle)




    def get_img(self):
        """returns the repartition of all pops"""

        out = np.ones((*self.dom.I.shape,3))

        out[:,:,0] = self.dom.I
        out[:,:,1] = self.dom.I
        out[:,:,2] = self.dom.I

        for s in self.states.values():
            out[np.where(self.Idx==s.idx)] = s.color


        bound = compute_bound(self.Idx)[0]
        out[np.where(bound == 1)] = np.array([0,0,0])


        return (out*255).astype(np.uint8)
        #return self.satisfaction<-12
        #return self.Idx
