import numpy as np

class Civ:
    def __init__(self, alpha, beta, d=2):
        self.alpha = alpha
        self.beta = beta
        self.d = d

    def c(self, P, R):
        return self.alpha*P/R

    def p(self, P, R):
        if self.alpha*P == 0:
            return 0
        return self.beta*R/P


class Parcel:
    def __init__(self, R0, Rmax, Pmin, Pmax, r):
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.R0 = R0
        self.Rmax = Rmax

        self.P = 0
        self.R = R0

        self.r = r
        self.civ = Civ(0,0,0)
        self.is_occupied = False
        self.is_overloaded = False
        self.is_collapsing = False
        self.is_occupable = True

    def block(self):
        self.is_occupable = False

    def set_civ(self, civ, Pinit):
        self.civ = civ
        self.P = Pinit
        self.is_occupied = True
        self.is_collapsing = False
        self.is_overloaded = False

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
        else:
            self.is_overloaded = False

    def add_resource(self, Radd):
        self.R += Radd
        if self.R >= self.Rmax:
            self.R = self.Rmax
        if self.R < 0:
            self.collapse()

    def update(self):
        if self.is_occupied:

            t_p = self.civ.p(self.P, self.R)
            t_c = self.r - self.civ.c(self.P, self.R)

            self.add_pop(t_p*self.P)
            self.add_resource(t_c*self.R)

            if t_p < 0 or t_c < 0:
                self.is_collapsing = True
            else:
                self.is_collapsing = False



class Simulation:
    def __init__(self, Map, Map_r, Map_R0, Map_Rmax, civ, pos_init=[2, 4]):
        self.Nx = Map.shape[0]
        self.Ny = Map.shape[1]
        self.start = np.array(pos_init)


        self.civ = civ

        self.parcels = []

        for i in range(self.Nx):
            self.parcels.append([])
            for j in range(self.Ny):

                if Map[i,j] == 0:
                    self.parcels[i].append(Parcel(0, 0, 0, 0, 0))
                    self.parcels[i][j].block()
                    self.parcels[i][j].P = -10
                else:
                    self.parcels[i].append(Parcel(Map_R0[i,j], Map_Rmax[i,j], 2, 100, Map_r[i,j]))
                    if [i,j] == pos_init:
                        self.parcels[i][j].set_civ(self.civ, 50)

    def occupied_cells(self):
        out = []
        for i in range(self.Nx):
            for j in range(self.Ny):
                if self.parcels[i][j].is_occupied:
                    out.append(np.array([i,j]))
        return np.array(out)


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
        bary = self.start
        ind = []
        for k in range(-self.civ.d, self.civ.d+1):
            for l in range(-self.civ.d, self.civ.d+1):
                if [k,l] != [0,0] and 0<=pos[0]+k<self.Nx and 0<=pos[1]+l<self.Ny:
                    if not(self.parcels[pos[0]+k][pos[1]+l].is_collapsing) and not(self.parcels[pos[0]+k][pos[1]+l].is_overloaded) and self.parcels[pos[0]+k][pos[1]+l].is_occupable:
                        dist = np.linalg.norm(np.array([pos[0]+k,pos[1]+l])-bary)
                        neighbors.append(self.parcels[pos[0]+k][pos[1]+l].R0/(dist**2))
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
                    best_pos = self.best_neighbor([i,j])
                    if best_pos:
                        trans_pop.append([[i,j], best_pos])
                elif self.parcels[i][j].is_collapsing:
                    best_pos = self.best_neighbor([i,j])
                    if best_pos:
                        if self.parcels[best_pos[0]][best_pos[1]].is_occupied:
                            trans_res.append([[i,j], best_pos])
                        else:
                            trans_pop.append([[i,j], best_pos])

        for ind in trans_res:
            self.transfer_res(*ind)

        for ind in trans_pop:
            self.transfer_pop(*ind)
