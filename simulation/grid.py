import numpy as np
from simulation.mathutils import grad_grad, lap


class Grid:
    def __init__(self, rho_0, pi_0, alpha, w, c, gamma, a, Kpi, r, Krho, DN_0, DR_0, bound):
        self.rho = rho_0
        self.pi = pi_0
        self.alpha = alpha
        self.w = w
        self.c = c
        self.gamma = gamma
        self.a = a
        self.Kpi = Kpi
        self.r = r
        self.Krho = Krho
        self.DN_0 = DN_0
        self.DR_0 = DR_0
        self.bound = bound
        self.area = np.sum(bound)

        self.repro = np.zeros_like(self.pi)

    def update(self):

        def DN(r, i):
            return self.DN_0[i]*np.exp(-r)

        for i in range(self.pi.shape[0]):
            conso = self.c[i]*np.sum([self.w[j]*self.a[i,j]*self.rho[j] for j in range(self.rho.shape[0])], axis=0)
            death = -self.gamma[i]
            war = -np.sum([self.alpha[i,k]*self.pi[k] for k in range(self.pi.shape[0])], axis=0)/self.Kpi[i]

            self.repro[i] = conso+death+war

            migration =  grad_grad(DN(self.repro[i], i),self.pi[i]) + DN(self.repro[i], i)*lap(self.pi[i])
            #migration = DN(self.pi[i], self.Kpi[i])*lap(self.pi[i])

            self.pi[i] += self.pi[i]*self.repro[i] + migration
            self.pi[i] *= self.bound

        def DR(r):
            return self.DR_0[i]*np.exp(r)

        for j in range(self.rho.shape[0]):
            renew = self.r[j]
            thresh = -(self.r[j]*self.rho[j])/self.Krho[j]
            conso = -np.sum([self.a[k,j]*self.pi[k] for k in range(self.pi.shape[0])], axis=0)

            repro_res = (renew+thresh+conso)
            exchange = np.sum([grad_grad(DR(self.repro[i]), self.rho[j]) + DR(self.repro[i])*lap(self.rho[j]) for i in range(self.pi.shape[0])], axis=0)
            #exchange = np.sum(DR(self), self.rho[j])*lap(self.rho[j]), axis=0)
            #exchange = 0

            self.rho[j] += self.rho[j]*repro_res + exchange
            self.rho[j] *= self.bound
