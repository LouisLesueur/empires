import numpy as np
from simulation.civ import Civ

class Parcel:
    def __init__(self, R0, Rmax, Pmin, Pmax, r, blanck_civ):
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.R0 = R0
        self.Rmax = Rmax
        self.blanck_civ = blanck_civ

        self.P = 0
        self.R = R0

        self.r = r
        self.civ = blanck_civ
        self.is_occupied = False
        self.is_overloaded = False
        self.is_collapsing = False
        self.is_at_war = False
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
        self.civ = self.blanck_civ
        self.is_occupied = False
        self.is_overloaded = False
        self.is_collapsing = False
        self.is_at_war = False
        self.is_occupable = True
        self.P = 0

    def add_pop(self, Padd):
        self.P += Padd
        if self.P >= self.Pmax:
            self.P = self.Pmax
            self.is_overloaded = True
        elif self.P <= self.Pmin:
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

            t_p = self.civ.p(self.P, self.R) - self.civ.gamma
            t_c = self.r - self.civ.c(self.P, self.R)

            if self.is_at_war:
                t_p -= 0.5
                t_c -= 0.5

            self.add_pop(t_p*self.P)
            self.add_resource(t_c*self.R)

            if t_p < 0 or t_c < 0:
                self.is_collapsing = True
            else:
                self.is_collapsing = False
