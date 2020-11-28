"""Grid module

This module combines domain, pops and resources to make the simulation
"""


import numpy as np
from grid.mathutils import lap, compute_bound, bound_lap


class Grid:
    """The grid class"""

    def __init__(self, N_start, start_pos, Ns, Rs, Domain, alpha, expend, dt):
        """Constructor of the grid

        N_start -- Starting pops
        start_pos -- Starting positions
        Ns -- pops
        Rs -- resources
        Domain -- The domain on which the simulation will be made
        alpha -- Initial values for each state
        expend -- Expension range (km)
        dt -- Time step (year)
        """

        self.dom = Domain
        self.dx = Domain.dx
        self.dt = dt

        self.Ns = Ns
        self.Rs = Rs

        self.r0 = Rs.r0*self.dom.I_topo*self.dom.I
        self.Rmax = Rs.Rmax*self.dom.I_topo
        self.Nmax = Ns.Nmax*self.dom.I_topo

        self.alpha = alpha
        self.expend = int(expend/self.dx)
        self.N = np.array(N_start, dtype=np.float64)

        self.citiesIdx = np.zeros_like(self.dom.I)-1
        self.Idx = np.zeros_like(self.dom.I)-1
        self.canExplore = np.zeros_like(self.dom.I, dtype=np.bool)

        self.city_pos = np.array(start_pos)

        self.neig = []
        for i in range(len(self.N)):
            self.neig.append([i])

            # id = [np.maximum(0,self.sites_pos[i,0]-self.expend),np.minimum(self.dom.I.shape[0],self.sites_pos[i,0]+self.expend),
            #       np.maximum(0,self.sites_pos[i,1]-self.expend),np.minimum(self.dom.I.shape[1],self.sites_pos[i,1]+self.expend)]
            #
            # self.citiesIdx[id[0]:id[1], id[2]:id[3]][self.Idx[id[0]:id[1], id[2]:id[3]]==-1] = i
            # self.Idx[id[0]:id[1], id[2]:id[3]][self.Idx[id[0]:id[1], id[2]:id[3]]==-1] = i

            self.citiesIdx[self.city_pos[i,0], self.city_pos[i,1]] = i
            self.canExplore[self.city_pos[i,0], self.city_pos[i,1]] = True
            self.Idx[self.city_pos[i,0], self.city_pos[i,1]] = i

        self.Idx[self.dom.I == 0] = -3

        self.R = self.Rmax[self.city_pos.T[0], self.city_pos.T[1]]
        self.Rpub = np.zeros_like(self.N)
        self.states = np.arange(len(self.N))
        self.colors = np.random.rand(len(self.states),3)

        self.time = 0

    def update(self):
        """Update the grid by one time step"""

        self.time += self.dt

        # Consumption
        # ---------------------------------------------------------------------
        self.conso = self.Ns.c0*(self.R/(self.Ns.Rdem + self.R))*self.N
        self.satisfaction = (self.R - self.Ns.Rdem)

        # Renewal
        # ---------------------------------------------------------------------
        renew = self.Rs.r0
        limit = -self.Rs.r0*(self.R/self.Rs.Rmax)

        dR = (renew+limit)*self.R-self.conso
        self.R += dR*self.dt

        # Taxes
        # ---------------------------------------------------------------------
        self.Rpub = self.alpha*self.R
        self.R = (1-self.alpha)*self.R

        self.satisfaction = (self.R - self.Ns.Rdem)

        # Demography (Bazykin model)
        # ---------------------------------------------------------------------

        death = - self.Ns.n0
        limit = - (self.N/self.Nmax[self.city_pos.T[0], self.city_pos.T[1]])
        G = death+limit

        dN = (G*self.N + self.Ns.chi * self.conso)*self.dt
        self.N += dN

        # Expension
        # ---------------------------------------------------------------------

        canstart = np.array(np.where(self.canExplore)).T
        start = canstart[np.random.randint(0,len(canstart), np.random.randint(0,len(canstart)))]
        oldlen = len(self.city_pos)

        for i, startpos in enumerate(start):
            city_idx = int(self.citiesIdx[startpos[0], startpos[1]])
            state_idx = int(self.Idx[startpos[0], startpos[1]])

            if np.random.rand() > 1-0.005:
                colony = startpos + np.random.randint(-1, 2, 2)*self.expend
                if (0 <colony[0] < self.dom.I.shape[0]) and (0 <colony[1] < self.dom.I.shape[1]):
                    if self.Idx[colony[0], colony[1]] == -1:
                        self.city_pos = np.append(self.city_pos, [colony], axis=0)
                        self.N = np.append(self.N, self.Ns.Nstart)
                        self.R = np.append(self.R, self.Rmax[colony[0], colony[1]])
                        self.neig[city_idx].append(oldlen)
                        self.neig.append([oldlen, city_idx])
                        self.citiesIdx[colony[0], colony[1]] = oldlen
                        self.Idx[colony[0], colony[1]] = state_idx
                        self.canExplore[colony[0], colony[1]] = True

            for id in range(-1,2):
                for jd in range(-1,2):
                    newpos = startpos + np.array([id,jd])
                    if (0 <newpos[0] < self.dom.I.shape[0]) and (0 <newpos[1] < self.dom.I.shape[1]):
                        if self.Idx[newpos[0], newpos[1]] == -1:
                            self.citiesIdx[newpos[0], newpos[1]] = city_idx
                            self.canExplore[newpos[0], newpos[1]] = True
                            self.Idx[newpos[0], newpos[1]] = state_idx
            self.canExplore[startpos[0], startpos[1]] = False


    def get_img(self):
        """returns the repartition of all pops"""

        out = np.ones((*self.dom.I.shape,3))

        out[:,:,0] = self.dom.I
        out[:,:,1] = self.dom.I
        out[:,:,2] = self.dom.I

        for i in range(len(self.states)):
            out[np.where(self.Idx==i)] = self.colors[i]


        bound = compute_bound(self.citiesIdx)[0]
        out[np.where(bound == 1)] = np.array([0,0,0])


        return (out*255).astype(np.uint8)
        #return self.satisfaction<-12
        #return self.canExplore
