"""Structure module

This module contains classes to represent political structures
"""

import matplotlib.pyplot as plt
import numpy as np
import json
from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image


class Domain:
    """The domain of integration"""

    def __init__(self, name):
        """Constructor of the domain

        name -- a string with the name of the domain files
        """

        f = open(name+'/specs.json',)
        specs = json.load(f)

        # Reading all the images
        self.I = plt.imread(name+'/bound.png')
        self.I_topo = plt.imread(name+'/topo.png')

        # Building the ressource base regeneration map
        self.I_r = (1-self.I_topo)*self.I

        self.dx = specs["dx"] #km
        self.area = np.sum(self.I)*(self.dx**2)
        self.shape = self.I.shape

    def dist(self, p1, p2):
        return np.linalg.norm(p1-p2)

    def exists(self, coords):
        if 0<=coords[0]<self.shape[0]:
            if 0<=coords[1]<self.shape[1]:
                return True
        return False




class States:

    def __init__(self, domain, max, max_reg, alpha, expend):

        self.number = 0

        self.dom = domain
        self.max = max
        self.max_reg = max_reg

        self.alpha = alpha

        self.res = np.zeros((max))-1
        self.relations = np.zeros((max,max))-1
        self.owned_regions = np.zeros((max_reg))-1

        self.expend = int(expend//self.dom.dx)

        self.colors = np.random.rand(max,3)

    #implement diplomacy and wars

    def taxes(self, region_res):
        for i,r in enumerate(self.owned_regions):
            if r>=0:
                self.res[int(r)] += self.alpha*region_res[i]
                region_res[i] -= self.alpha*region_res[i]

    def add_colony(self, id_from, id_to):
        state_from = self.owned_regions[id_from]
        self.owned_regions[id_to] = state_from

    def add_state(self, init_region):

        if self.number < self.max:
            self.owned_regions[init_region] = self.number
            self.number += 1

    def update_diplomacy(self, roads, margin):
        # 0: pacific
        # 1: neutral
        # 2: agressive
        for i in range(self.max):
            for j in range(self.max):
                if np.abs(self.res[i] - self.res[j]) < margin:
                    self.relations[i,j] = 1
                    self.relations[j,i] = 1
                elif self.res[i] < self.res[j]:
                    self.relations[i,j] = 0
                    self.relations[j,i] = 2
                else:
                    self.relations[i,j] = 2
                    self.relations[j,i] = 0





class Regions:
    def __init__(self, states):

        self.number = 0

        self.states = states
        self.max = self.states.max_reg

        self.map = self.states.dom.I.copy() - 2
        self.canExplore = np.zeros_like(self.states.dom.I, dtype=np.bool)

        self.roads = np.zeros((self.max, self.max)) #1: road, 2: contact
        self.cities = np.zeros((self.max,2))
        self.pop = np.zeros(self.max)
        self.res = np.zeros(self.max)
        self.r = np.zeros(self.max)
        self.area = np.zeros(self.max)
        self.stab = 1000+np.zeros(self.max)


    def revolts(self):
        if np.random.rand() > 1-0.05:
            for city in np.where(self.stab < 0)[0]:
                if np.random.rand() > 1-0.05:
                    self.states.add_state(city)

    def wars_and_alliances(self):
        if np.random.rand() > 1-0.1:
            cancoop = np.array(np.where(self.roads == 2)).T

            if len(cancoop) > 0:
                coop = cancoop[np.random.randint(0,len(cancoop), np.random.randint(0,len(cancoop)))]


                for pair in coop:
                    if self.states.relations[pair[0], pair[1]] == 0:
                        self.roads[pair[0], pair[1]] = 1
                        self.roads[pair[1], pair[0]] = 1
                    elif self.states.relations[pair[0], pair[1]] == 2:
                        self.states.owned_regions[pair[1]] = self.states.owned_regions[pair[0]]
                        self.roads[pair[0], pair[1]] = 1
                        self.roads[pair[1], pair[0]] = 1




    def update(self, c0, chi, n0, Rm, dt):
        resMax = self.area*10 + 50

        popMax = self.area*100 + 1e-6

        conso = c0*(self.res/(Rm + self.res))*self.pop
        self.stab[self.pop>0] = self.res[self.pop>0] - Rm

        renew = self.r
        limit = -self.r*(self.res/resMax)

        dres = (renew+limit)*self.res
        self.res += (dres-conso)*dt

        self.states.taxes(self.res)

        death = -n0
        limit = -self.pop/popMax
        G = death+limit

        dN = (G*self.pop + chi * conso)*dt
        self.pop += dN

        self.pop = np.maximum(0,self.pop)
        self.res = np.maximum(0,self.res)



    def add_city_state(self, pos, init_pop, init_res):

        if self.states.dom.exists(pos):
            if self.map[pos[0], pos[1]] == -1:
                self.cities[self.number] = pos
                self.pop[self.number] = init_pop
                self.res[self.number] = init_res
                self.r[self.number] = self.states.dom.I_r[pos[0], pos[1]]
                self.area[self.number] = 1
                self.canExplore[pos[0], pos[1]] = True
                self.map[pos[0], pos[1]] = self.number
                self.states.add_state(self.number)
                self.number += 1
                return True
        return False


    def colonize(self, coords_from, id_from, coords_to):

        if self.number > self.max:
            return

        id_to = self.number

        #coords_from = self.cities[id_from]

        dist = int(self.states.dom.dist(coords_from, coords_to))
        path = np.array([coords_from*(i/dist)+coords_to*(1-(i/dist)) for i in range(dist)], dtype=np.int32)
        pathidx = np.array([self.map[rp[0], rp[1]] for rp in path])

        if not(len(pathidx[pathidx > 0])>0):
            for i,p in enumerate(path):
                if self.map[p[0], p[1]] == -1:
                    if i >= len(path)//2:
                        self.map[p[0], p[1]] = id_from
                        self.r[id_from] += self.states.dom.I_r[p[0], p[1]]
                        self.area[id_from] += 1
                    else:
                        self.map[p[0], p[1]] = id_to
                        self.r[id_to] += self.states.dom.I_r[p[0], p[1]]
                        self.area[id_to] += 1
                    self.canExplore[p[0], p[1]] = True

            self.roads[id_from, id_to] = 1
            self.roads[id_to, id_from] = 1

            self.cities[id_to] = coords_to
            self.pop[id_to] = 0.2*self.pop[id_from]
            self.res[id_to] = 0.2*self.res[id_from]
            self.r[id_to] = self.states.dom.I_r[coords_to[0], coords_to[1]]
            self.states.add_colony(id_from, id_to)

            self.number += 1


    def grow(self):

        canstart = np.array(np.where(self.canExplore)).T
        start = canstart[np.random.randint(0,len(canstart), np.random.randint(0,len(canstart)))]

        for i, startpos in enumerate(start):
            city_idx = int(self.map[startpos[0], startpos[1]])

            if np.random.rand() > 1-0.03:
                go_to = startpos + np.random.randint(-1, 2, 2)*self.states.expend
                if self.states.dom.exists(go_to):
                    if self.map[go_to[0], go_to[1]] == -1:
                        self.colonize(startpos, city_idx, go_to)

            for id in range(-1,2):
                for jd in range(-1,2):
                    newpos = startpos + np.array([id,jd])
                    if self.states.dom.exists(newpos):
                        if int(self.map[newpos[0], newpos[1]]) == -1:
                            self.map[newpos[0], newpos[1]] = city_idx
                            self.r[city_idx] += self.states.dom.I_r[newpos[0], newpos[1]]
                            self.area[self.number] += 1
                            self.canExplore[newpos[0], newpos[1]] = True
                        elif int(self.map[newpos[0], newpos[1]]) > 0 and self.states.owned_regions[int(self.map[newpos[0], newpos[1]])] != self.states.owned_regions[city_idx]:
                            self.roads[int(self.map[newpos[0], newpos[1]]), city_idx] = 2
                            self.roads[city_idx, int(self.map[newpos[0], newpos[1]])] = 2

            self.canExplore[startpos[0], startpos[1]] = False
