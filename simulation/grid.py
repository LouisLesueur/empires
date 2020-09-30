import numpy as np
from simulation.mathutils import grad_grad, lap


class Grid:
    def __init__(self, rho_0, pi_0, alpha,  gamma, a, Kpi, r, K, DN_0, bound):
        self.rho = rho_0
        self.pi = pi_0
        self.alpha = alpha
        self.gamma = gamma

        self.a = a
        self.r = r
        self.K = K
        self.Kpi = Kpi
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
            war = 0

            self.repro[i] = conso+death+war

            migration =  DN(self.repro[i], i)*lap(self.pi[i])
            shift = 0.001*self.pi[i]*lap(np.sum(self.rho, axis=0))

            self.pi[i] += self.pi[i]*self.repro[i] + migration+shift
            self.pi[i] *= self.bound


        for j in range(self.rho.shape[0]):
            renew = self.r[j]
            thresh = -(self.r[j]*self.rho[j])/self.K[j]
            conso = -np.sum([self.a[k,j]*self.pi[k] for k in range(self.pi.shape[0])], axis=0)

            repro_res = (fluctuations*renew+thresh+conso)

            self.rho[j] += self.rho[j]*repro_res
            self.rho[j] *= self.bound


        def masking(i):
            mask = np.zeros_like(self.pi[i])
            mask[np.where(self.pi[i]>0.01)]=1
            return mask

        for i in range(self.pi.shape[0]):
            for k in range(self.pi.shape[0]):
                if ( i != k):
                    PI = [np.sum(self.pi[k]) for k in range(self.pi.shape[0])]
                    self.alpha[i,k] = PI[k]/self.Kpi[i]
            self.alpha[i,:] /= np.linalg.norm(self.alpha[i,:])


            for j in range(self.rho.shape[0]):
                lambada = 0.00001
                PI = [np.sum(self.pi[k])*self.area for k in range(self.pi.shape[0])]
                R = [np.sum(self.r[k])*self.area for k in range(self.rho.shape[0])]
                K_mean = np.sum(self.K[j]*self.area)

                self.a[i,j] = R[j]/(2*lambada*K_mean)
