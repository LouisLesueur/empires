import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


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
    def __init__(self,R0, N0, alpha, w, c, gamma, a, KN, r, KR, bound):
        self.R = R0
        self.N = N0
        self.alpha = alpha
        self.w = w
        self.c = c
        self.gamma = gamma
        self.a = alpha
        self.KN = KN
        self.r = r
        self.KR = KR
        self.bound = bound

    def update(self):

        def DN(N, K):
            out = np.sum([Rr for Rr in self.R], axis=0)
            out[np.where(2*N/K > 1)] = 0
            for i in range(self.R.shape[0]):
                out[np.where(2*self.R[i]/self.KR[i]<1)]=0
            return out

        for i in range(self.N.shape[0]):
            conso = self.c[i]*np.sum([self.w[j]*self.a[i,j]*self.R[j] for j in range(self.R.shape[0])], axis=0)
            death = -self.gamma[i]
            war = - np.sum([self.alpha[i,k]*self.N[k] for k in range(self.N.shape[0])], axis=0)/self.KN[i]

            migration =  grad_grad(DN(self.N[i], self.KN[i]), self.N[i]) + DN(self.N[i], self.KN[i])*lap(self.N[i])
            #migration = 0.00005*DN(self.N[i], self.KN[i])*lap(self.N[i])

            self.N[i] += self.N[i]*(conso+death+war) + migration
            self.N[i] *= self.bound

        def DR(R, K):
            out = np.sum([Nr for Nr in self.N], axis=0)
            out[np.where(2*R/K < 1)] = 0
            for i in range(self.N.shape[0]):
                out[np.where(2*self.N[i]/self.KN[i]>1)]=0
            return out

        for j in range(self.R.shape[0]):
            renew = self.r[j]
            thresh = -(self.r[j]*self.R[j])/self.KR[j]
            conso = -np.sum([self.a[k,j]*self.N[k] for k in range(self.N.shape[0])], axis=0)

            exchange = grad_grad(DR(self.R[j], self.KR[j]), self.R[j]) + DR(self.R[j], self.KR[j])*lap(self.R[j])
            #exchange = 0.00005*DR(self.R[j], self.KR[j])*lap(self.R[j])

            self.R[j] += self.R[j]*(renew+thresh+conso) + exchange
            self.R[j] *= self.bound



I = plt.imread('maps/europe.png')
I_r = plt.imread('maps/europe_r.png')
I_R0 = plt.imread('maps/europe_r.png')

I = 0.2989 * I[:,:,0] + 0.5870 * I[:,:,1] + 0.1140 * I[:,:,2]
I_r = (0.2989 * I_r[:,:,0] + 0.5870 * I_r[:,:,1] + 0.1140 * I_r[:,:,2])*0.08
I_R0 = (0.2989 * I_R0[:,:,0] + 0.5870 * I_R0[:,:,1] + 0.1140 * I_R0[:,:,2])*1250+1

N1 = np.zeros_like(I)
N2 = np.zeros_like(I)

N1[50, 62] = 20  #Rome
N2[20, 120] = 20 #Barbarians


sim = Grid(np.array([I_R0]), #R
           np.array([N1,N2]),
           np.array([[0,0.04],  #alpha
                     [0.03,0]]),
           np.array([1]), #w
           np.array([0.01, 0.01]), #c
           np.array([0.05*np.ones(N1.shape), 0.05*np.ones(N2.shape)]), #gamma
           np.array([[1],
                     [1]]), #a
           np.array([1000, 400]), #KN
           np.array([I_r]), #r
           np.array([I_R0]),
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

    N_out = np.zeros((N.shape[0], N.shape[1], N.shape[2], 3))
    for i in range(len(N)):
        N_out[i] = I
        N_out[i][np.where(N[i]>0)] = colors[0]


fig = plt.figure()
ax = plt.axes()
im=plt.imshow(sim.N[0]+sim.N[1],interpolation='none')

def animate(i):
    sim.update()
    im.set_array(sim.N[0]+sim.N[1])
    return [im]

writer = PillowWriter(fps=25)
anim = FuncAnimation(fig, animate, frames = 500, interval = 100)
#anim.save('out.gif', writer=writer)
plt.colorbar()
plt.show()
