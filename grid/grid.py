"""Grid module

This module combines domain, pops and resources to make the simulation
"""


import numpy as np
from grid.mathutils import lap, compute_bound, bound_lap

class Simulation:
    def __init__(self, regions, c0, chi, n0, Rm, dt):
        self.regions = regions
        self.dom = self.regions.states.dom
        self.c0 = c0
        self.chi=chi
        self.n0 = n0
        self.Rm = Rm
        self.dt = dt

        self.time = 0

    def random_city_state(self, N):

        count = 0
        while count<N:
            i = np.random.randint(10, self.dom.shape[0]-10)
            j = np.random.randint(10, self.dom.shape[1]-10)

            pos = np.array([i,j])
            if self.regions.add_city_state(pos, 5, 5):
                count += 1

    def update(self):
        self.time += self.dt
        self.regions.update(self.c0, self.chi, self.n0, self.Rm, self.dt)
        self.regions.grow()

        if self.time > 5*self.dt:
            #self.regions.revolts()
            self.regions.states.update_diplomacy(self.regions.roads, 0.1)
            self.regions.wars_and_alliances()

    def get_img(self):
        out = np.ones((*self.dom.I.shape,3))

        out[:,:,0] = self.dom.I
        out[:,:,1] = self.dom.I
        out[:,:,2] = self.dom.I

        for i,reg in enumerate(self.regions.states.owned_regions):
            out[np.where(self.regions.map==i)] = self.regions.states.colors[int(reg)]


        return (out*255).astype(np.uint8)

    def get_roads(self):
        roads = []
        potential = []
        for i in range(self.regions.max):
            for j in range(i,self.regions.max):
                if self.regions.roads[i,j] == 1:
                    roads.append([(self.regions.cities[i][1], self.regions.cities[i][0]),
                                  (self.regions.cities[j][1], self.regions.cities[j][0])])
                elif self.regions.roads[i,j] == 2:
                    potential.append([(self.regions.cities[i][1], self.regions.cities[i][0]),
                                      (self.regions.cities[j][1], self.regions.cities[j][0])])
        return roads, potential
