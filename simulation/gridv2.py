import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def grad(D, axis=0):
    D_pad = np.pad(D, (1,1), 'reflect')
    return (np.roll(D_pad,1,axis=axis) - D_pad)[1:-1,1:-1]

def grad_grad(D, A):
    A_pad = np.pad(A, (1,1), 'reflect')
    D_pad = np.pad(D, (1,1), 'reflect')

    D_x = (np.roll(D_pad,1,axis=0) - D_pad)[1:-1,1:-1]
    A_x = (np.roll(A_pad,1,axis=0) - A_pad)[1:-1,1:-1]

    D_y = (np.roll(D_pad,1,axis=1) - D_pad)[1:-1,1:-1]
    A_y = (np.roll(A_pad,1,axis=1) - A_pad)[1:-1,1:-1]

    return D_x*A_x + D_y*A_y

def lap(A):
    A_pad = np.pad(A, (1,1), 'reflect')
    return ((np.roll(A_pad,-1,axis=1) - 2*A_pad + np.roll(A_pad,1,axis=1)) + (np.roll(A_pad,-1,axis=0) - 2*A_pad + np.roll(A_pad,1,axis=0)))[1:-1,1:-1]



class Grid:
    def __init__(self, rho_0, pi_0, alpha, w, c, gamma, a, Kpi, r, Krho, bound):
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
        self.bound = bound
        self.area = np.sum(bound)

        self.repro = np.zeros_like(self.pi)

    def update(self):

        def DN(r, rr):
            res = np.sum(ro for ro in self.rho)
            return (res/np.linalg.norm(res))*np.exp(-r)

        for i in range(self.pi.shape[0]):
            conso = self.c[i]*np.sum([self.w[j]*self.a[i,j]*self.rho[j] for j in range(self.rho.shape[0])], axis=0)
            death = -self.gamma[i]
            war = -np.sum([self.alpha[i,k]*self.pi[k] for k in range(self.pi.shape[0])], axis=0)/self.Kpi[i]

            self.repro[i] = conso+death+war

            migration =  grad_grad(DN(self.repro[i], np.sum(ro for ro in self.rho)),self.pi[i]) + DN(self.repro[i], np.sum(ro for ro in self.rho))*lap(self.pi[i])
            #migration = DN(self.pi[i], self.Kpi[i])*lap(self.pi[i])

            self.pi[i] += self.pi[i]*self.repro[i] + migration
            self.pi[i] *= self.bound

        def DR(r):
            return 0.02*np.exp(0.01*r)

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



I = plt.imread('maps/europe.png')
I_r = plt.imread('maps/europe_r.png')
I_R0 = plt.imread('maps/europe_r.png')

I = 0.2989 * I[:,:,0] + 0.5870 * I[:,:,1] + 0.1140 * I[:,:,2]
I_r = (0.2989 * I_r[:,:,0] + 0.5870 * I_r[:,:,1] + 0.1140 * I_r[:,:,2])*0.08
I_R0 = (0.2989 * I_R0[:,:,0] + 0.5870 * I_R0[:,:,1] + 0.1140 * I_R0[:,:,2])*100+1

N1 = np.zeros_like(I)
N2 = np.zeros_like(I)

N1[50, 62] = 20  #Rome
N2[20, 120] = 20 #Barbarians


sim = Grid(np.array([I_R0]), #R
           np.array([N1,N2]),
           np.array([[0,0.4],  #alpha
                     [0.03,0]]),
           np.array([1]), #w
           np.array([0.2, 0.1]), #c
           np.array([0.01*np.ones(N1.shape), 0.05*np.ones(N2.shape)]), #gamma
           np.array([[0.005],
                     [0.01]]), #a
           np.array([30, 40]), #KN
           np.array([I_r]), #r
           np.array([2*I_R0]), #KR
           I) #KR

#Tmax = 1500
#N_tot1 = np.zeros(Tmax)
#N_tot2 = np.zeros(Tmax)
#R_tot = np.zeros(Tmax)

#N_tot1[0]=np.sum(sim.N[0])
#N_tot2[0]=np.sum(sim.N[1])
#R_tot[0]=np.sum(sim.N[1])

#for i in range(Tmax):
#    sim.update()
#    N_tot1[i]=np.sum(sim.N[0])
#    N_tot2[i]=np.sum(sim.N[1])
#    R_tot[i]=np.sum(sim.R[0])

#plt.plot(R_tot, label='Ressources')


def colormap(I, N, colors):

    N_out = I.copy()
    for i in range(len(N)):
        for k in range(N[i].shape[0]):
            for l in range(N[i].shape[1]):
                if N[i][k,l] > 0:
                    N_out[k,l] += colors[i]
    return N_out[:,:,:3]


I = plt.imread('maps/europe.png')
print(I)


good = np.zeros_like(sim.pi[0])
good[np.where(sim.pi[0]>5)] = 1
bad = np.zeros_like(sim.pi[0])
bad[np.where(sim.pi[1]>5)] = 2
fig = plt.figure()
ax = plt.axes()
im=plt.imshow(good+bad)

def animate(i):

    good = np.zeros_like(sim.pi[0])
    good[np.where(sim.pi[0]>5)] = 1
    bad = np.zeros_like(sim.pi[0])
    bad[np.where(sim.pi[1]>5)] = 2
    sim.update()
    im.set_array(good+bad)
    return [im]

writer = PillowWriter(fps=25)
anim = FuncAnimation(fig, animate, frames = 500, interval = 100)
#anim.save('out.gif', writer=writer)
plt.colorbar()
plt.show()
