import numpy as np

class Civ:
    def __init__(self, alpha, beta, d=2):
        self.alpha = alpha
        self.beta = beta
        self.d = d

    def c(self, P, Pmax):
        return self.alpha*P/Pmax

    def p(self, R, P):
        if self.alpha*P == 0:
            return 0
        return self.beta*np.tanh(R/(self.alpha*P)-1)


class Parcel:
    def __init__(self, R0, Pmin, Pmax, r):
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.R0 = R0

        self.P = 0
        self.R = R0

        self.r = r
        self.civ = Civ(0,0,0)
        self.is_occupied = False
        self.is_overloaded = False
        self.is_collapsing = False

    def set_civ(self, civ, Pinit):
        self.civ = civ
        self.P = Pinit
        self.is_occupied = True

    def collapse(self):
        self.civ = Civ(0,0,0)
        self.is_occupied = False
        self.P = 0

    def add_pop(self, Padd):
        self.P += Padd
        if self.P >= self.Pmax:
            self.P = self.Pmax
            self.is_overloaded = True
        elif self.P < self.Pmin:
            self.collapse()

    def add_resource(self, Radd):
        self.R += Radd
        if self.R < 0:
            self.collapse()

    def update(self):
        if self.is_occupied:
            self.add_pop(self.civ.p(self.R, self.P)*self.P)
            self.add_resource((self.r - self.civ.c(self.P, self.Pmax))*self.R)
            if self.civ.p(self.R, self.P) < 0:
                self.is_collapsing = True



class Grid:
    def __init__(self, Nx, Ny, pos_init=[0,0]):

        self.Nx = Nx
        self.Ny = Ny

        self.pos_init = pos_init

        self.civ = Civ(0.2, 0.2, 2)

        self.parcels = []

        for i in range(Nx):
            self.parcels.append([])
            for j in range(Ny):
                self.parcels[i].append(Parcel(np.random.randint(10,30), 20, 200, np.random.random()))
                if [i,j] == pos_init:
                    self.parcels[i][j].set_civ(self.civ, 50)

    def transfer_pop(self, pos_from, pos_to, fact=0.5):
        Ptot = fact*self.parcels[pos_from[0]][pos_from[1]].P
        self.parcels[pos_from[0]][pos_from[1]].add_pop(-Ptot)
        if self.parcels[pos_to[0]][pos_to[1]].is_occupied:
            self.parcels[pos_to[0]][pos_to[1]].add_pop(Ptot)
        else:
            self.parcels[pos_to[0]][pos_to[1]].set_civ(self.parcels[pos_from[0]][pos_from[1]].civ, Ptot)

    def transfer_res(self, pos_from, pos_to, fact=0.25):
        Rtot = fact*self.parcels[pos_from[0]][pos_from[1]].R
        self.parcels[pos_from[0]][pos_from[1]].add_resource(-Rtot)
        self.parcels[pos_to[0]][pos_to[1]].add_resource(Rtot)


    def best_neighbor(self, pos):
        neighbors = []
        ind = []
        for k in range(-self.civ.d, self.civ.d+1):
            for l in range(-self.civ.d, self.civ.d+1):
                if [k,l] != [0,0] and 0<=pos[0]+k<self.Nx and 0<=pos[1]+l<self.Ny:
                    if not(self.parcels[pos[0]+k][pos[1]+l].is_collapsing):
                        neighbors.append(self.parcels[pos[0]+k][pos[1]+l].R0)
                        ind.append([pos[0]+k,pos[1]+l])
        if neighbors != []:
            idx = np.argmax(neighbors)
            return ind[idx]

        else:
            return False


    def update(self):
        for i in range(self.Nx):
            for j in range(self.Ny):
                self.parcels[i][j].update()
                if self.parcels[i][j].is_overloaded:
                    best_pos = self.best_neighbor([i,j])
                    if best_pos:
                        self.transfer_pop([i,j], best_pos)
                if self.parcels[i][j].is_collapsing:
                    best_pos = self.best_neighbor([i,j])
                    if best_pos:
                        if self.parcels[best_pos[0]][best_pos[1]].is_occupied:
                            self.transfer_res([i,j], best_pos)
