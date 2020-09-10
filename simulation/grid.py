import numpy as np
from simulation.civ import Civ
from simulation.parcel import Parcel

class Grid:
    def __init__(self, Map, Map_r, Map_R0, Map_Rmax, civ):
        self.Nx = Map.shape[0]
        self.Ny = Map.shape[1]

        greenCiv = Civ(0, 0, 0, 0 ,None, np.array([0,1,0]))
        waterCiv = Civ(0, 0, 0, 0 ,None, np.array([0,1,1]))

        blanckParcel = Parcel(0, 0, 0, 0, 0, waterCiv)
        blanckParcel.block()
        blanckParcel.P = -100

        self.parcels = []

        for i in range(self.Nx):
            self.parcels.append([])
            for j in range(self.Ny):

                if Map[i,j] == 0:
                    self.parcels[i].append(blanckParcel)
                else:
                    self.parcels[i].append(Parcel(Map_R0[i,j], Map_Rmax[i,j], 1, 10000, Map_r[i,j], greenCiv))
                    for k in range(len(civ)):
                        if [i,j] == civ[k].capital:
                            self.parcels[i][j].set_civ(civ[k], 50)


    def transfer_pop(self, pos_from, pos_to, fact=0.5):
        Ptot = fact*self.parcels[pos_from[0]][pos_from[1]].P

        if self.parcels[pos_to[0]][pos_to[1]].is_occupied:
            if self.parcels[pos_to[0]][pos_to[1]].civ == self.parcels[pos_from[0]][pos_from[1]].civ:
                self.parcels[pos_from[0]][pos_from[1]].add_pop(-Ptot)
                self.parcels[pos_to[0]][pos_to[1]].add_pop(Ptot)
                self.parcels[pos_to[0]][pos_to[1]].is_at_war = False
                self.parcels[pos_from[0]][pos_from[1]].is_at_war = False
            else:
                self.parcels[pos_to[0]][pos_to[1]].is_at_war = True
                self.parcels[pos_from[0]][pos_from[1]].is_at_war = True
        else:
            self.parcels[pos_from[0]][pos_from[1]].add_pop(-Ptot)
            self.parcels[pos_to[0]][pos_to[1]].set_civ(self.parcels[pos_from[0]][pos_from[1]].civ, Ptot)

    def transfer_res(self, pos_from, pos_to, fact=0.25):
        Rtot = fact*self.parcels[pos_from[0]][pos_from[1]].R

        if self.parcels[pos_to[0]][pos_to[1]].civ == self.parcels[pos_from[0]][pos_from[1]].civ:
            self.parcels[pos_from[0]][pos_from[1]].add_resource(-Rtot)
            self.parcels[pos_to[0]][pos_to[1]].add_resource(Rtot)
            self.parcels[pos_to[0]][pos_to[1]].is_at_war = False
            self.parcels[pos_from[0]][pos_from[1]].is_at_war = False
        else:
            self.parcels[pos_to[0]][pos_to[1]].is_at_war = True
            self.parcels[pos_from[0]][pos_from[1]].is_at_war = True



    def best_neighbor(self, pos, d):
        neighbors = []
        bary = self.parcels[pos[0]][pos[1]].civ.capital
        ind = []
        for k in range(-d, d+1):
            for l in range(-d, d+1):
                if [k,l] != [0,0] and 0<=pos[0]+k<self.Nx and 0<=pos[1]+l<self.Ny:
                    if not(self.parcels[pos[0]+k][pos[1]+l].is_collapsing) and not(self.parcels[pos[0]+k][pos[1]+l].is_overloaded) and self.parcels[pos[0]+k][pos[1]+l].is_occupable:
                        dist = np.linalg.norm(np.array([pos[0]+k,pos[1]+l])-bary)
                        neighbors.append(self.parcels[pos[0]+k][pos[1]+l].R/(dist**2))
                        ind.append([pos[0]+k,pos[1]+l])
                    elif self.parcels[pos[0]+k][pos[1]+l].is_occupied and self.parcels[pos[0]+k][pos[1]+l].civ != self.parcels[pos[0]][pos[1]].civ:
                        dist = np.linalg.norm(np.array([pos[0]+k,pos[1]+l])-bary)
                        neighbors.append(self.parcels[pos[0]+k][pos[1]+l].R/(dist**2))
                        ind.append([pos[0]+k,pos[1]+l])

        if neighbors != []:
            idx = np.argmax(neighbors)
            return ind[idx]

        else:
            return False


    def update(self):
        trans_res = []
        trans_pop = []

        for i in range(self.Nx):
            for j in range(self.Ny):
                self.parcels[i][j].update()

                if self.parcels[i][j].is_overloaded:
                    best_pos = self.best_neighbor([i,j], self.parcels[i][j].civ.d)
                    if best_pos:
                        trans_pop.append([[i,j], best_pos])

                elif self.parcels[i][j].is_collapsing:
                    best_pos = self.best_neighbor([i,j], self.parcels[i][j].civ.d)
                    if best_pos:
                        if self.parcels[best_pos[0]][best_pos[1]].is_occupied:
                            trans_res.append([[i,j], best_pos])
                        else:
                            trans_pop.append([[i,j], best_pos])

        for ind in trans_res:
            self.transfer_res(*ind)

        for ind in trans_pop:
            self.transfer_pop(*ind)
