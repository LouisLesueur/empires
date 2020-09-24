import matplotlib.pyplot as plt
import numpy as np
from simulation.grid import Grid

from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image


I = plt.imread('maps/europe.png')
I_r = plt.imread('maps/europe_r.png')
I_R0 = plt.imread('maps/europe_r.png')

I = 0.2989 * I[:,:,0] + 0.5870 * I[:,:,1] + 0.1140 * I[:,:,2]
I_r = (0.2989 * I_r[:,:,0] + 0.5870 * I_r[:,:,1] + 0.1140 * I_r[:,:,2])*0.08
I_R0 = (0.2989 * I_R0[:,:,0] + 0.5870 * I_R0[:,:,1] + 0.1140 * I_R0[:,:,2])*100+1

N1 = np.zeros_like(I)
N2 = np.zeros_like(I)
N3 = np.zeros_like(I)

N1[50, 62] = 20  #Rome
N2[20, 120] = 20 #Barbarians
N3[34, 44] = 20  #Rome


sim = Grid(np.array([I_R0]), #R
           np.array([N1,N2,N3]),
           np.array([[0,0.4,0.2],  #alpha
                     [0.03,0,0.02],
                     [0,0.2,0.5]]),
           np.array([1]), #w
           np.array([0.15, 0.15, 0.15]), #c
           np.array([0.02*np.ones(N1.shape), 0.05*np.ones(N2.shape), 0.03*np.ones(N2.shape)]), #gamma
           np.array([[0.005],
                     [0.01],
                     [0.007]]), #a
           np.array([30, 40, 50]), #KN
           np.array([I_r]), #r
           np.array([2*I_R0]), #KR
           np.array([0.01, 0.01, 0.0005]), #DN_0
           np.array([2, 10, 0.005]), #DR_0
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


def colormap(N, thresh, color):
    good = np.zeros_like(N)
    good[np.where(N>thresh)] = color
    return good


I = plt.imread('maps/europe.png')

threshold = 4

#out = np.zeros_like(sim.pi[0])
#for i in range(len(sim.pi)):
#    out += colormap(sim.pi[i], threshold, i+1)
fig = plt.figure()
ax = plt.axes()
im=plt.imshow(np.sum(sim.pi, axis=0),cmap=plt.cm.magma)

def animate(i):

    #out = np.zeros_like(sim.pi[0])
    #for i in range(len(sim.pi)):
    #    out += colormap(sim.pi[i], threshold, i+1)
    sim.update()
    im.set_array(np.sum(sim.pi, axis=0))
    return [im]

writer = PillowWriter(fps=25)
anim = FuncAnimation(fig, animate, frames = 1000, interval = 100)
#anim.save('out.gif', writer=writer)
plt.colorbar()
plt.show()
