import numpy as np
from simulation.mathutils import grad_grad, lap


class Grid:
    def __init__(self, rho_0, pi_0,  gamma, a, r, K, DN_0, bound):
        self.rho = rho_0
        self.pi = pi_0
        self.gamma = gamma

        self.a = a
        self.r = r
        self.K = K
        self.DN_0 = DN_0
        self.bound = bound
        self.area = np.sum(bound)

        self.time = 0
        self.repro = np.zeros_like(self.pi)


    def update(self):

        self.time += 1
        fluctuations = 2*(1-np.sin(100*self.time))

        def DN(r, i):
            return self.DN_0[i]*np.exp(-(r/np.linalg.norm(r)))

        for i in range(self.pi.shape[0]):
            conso = np.sum([self.a[i,j]*self.rho[j] for j in range(self.rho.shape[0])], axis=0)
            death = -self.gamma[i]

            self.repro[i] = conso+death

            migration =  DN(self.repro[i], i)*lap(self.pi[i])
            shift = 0.001*self.pi[i]*lap(np.sum(self.rho, axis=0))

            self.pi[i] += self.pi[i]*self.repro[i] + migration
            self.pi[i] *= self.bound

        for j in range(self.rho.shape[0]):
            renew = self.r[j]
            thresh = -(self.r[j]*self.rho[j])/self.K[j]
            conso = -np.sum([self.a[k,j]*self.pi[k] for k in range(self.pi.shape[0])], axis=0)

            repro_res = (fluctuations*renew+thresh+conso)

            self.rho[j] += self.rho[j]*repro_res
            self.rho[j] *= self.bound



        for i in range(self.pi.shape[0]):
            for j in range(self.rho.shape[0]):
                lambada = 0.01
                fac = 0.5*(1/(lambada + self.K[j]*(self.pi[i]/(self.r[j]+1e-6))))
                self.a[i,j] = fac*(self.K[j]-np.sum([(self.K[j])*self.a[k,j]*self.pi[k]/(self.r[j]+1e-6) for k in range(self.rho.shape[0]) if k != i]))
                self.a[i,j] *= self.bound
                self.a[i,j] /= np.linalg.norm(self.a[i,j])
