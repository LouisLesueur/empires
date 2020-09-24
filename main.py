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
N4 = np.zeros_like(I)

N1[50, 62] = 20  #Rome
N2[30, 83] = 20 #Barbarians
N3[34, 44] = 20  #Paris
N4[25, 100] = 20  #Russe


sim = Grid(np.array([I_R0]), #R
           np.array([N1,N2,N3,N4]),
           np.array([[0,   0.4,0.2, 0.1],  #alpha
                     [0.03,0,  0.02, 0.01],
                     [0.04, 0.2, 0, 0.5],
                     [0.01, 0.04, 0.05, 0]]),
           np.array([1]), #w
           np.array([0.15, 0.15, 0.15, 0.15]), #c
           np.array([0.02*np.ones(N1.shape), 0.05*np.ones(N2.shape), 0.03*np.ones(N2.shape), 0.02*np.ones(N2.shape)]), #gamma
           np.array([[0.005],
                     [0.006],
                     [0.007],
                     [0.008]]), #a
           np.array([30, 40, 50, 40]), #KN
           np.array([I_r]), #r
           np.array([2*I_R0]), #KR
           np.array([0.01, 0.01, 0.0005, 0.004]), #DN_0
           np.array([2, 10, 0.005, 0.01]), #DR_0
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



def colormap(pi, color):
    out = color*np.ones((pi.shape[0], pi.shape[1], 3))
    fact = pi/np.max(pi)
    out[:,:,0] *= fact
    out[:,:,1] *= fact
    out[:,:,2] *= fact
    return out


I = plt.imread('maps/europe.png')

colors = np.array([[0,0,1],[0,1,0],[1,0,0], [1,0,1]])
out = np.zeros((sim.pi[0].shape[0], sim.pi[0].shape[1], 3))
for i in range(len(sim.pi)):
    out += colormap(sim.pi[i], colors[i])

fig = plt.figure()
ax = plt.axes()
im=plt.imshow(out)


def animate(i):

    out = np.zeros((sim.pi[0].shape[0], sim.pi[0].shape[1], 3))
    for k in range(len(sim.pi)):
        out += colormap(sim.pi[k], colors[k])
    for l in range(7):
        sim.update()
    im.set_array(out)
    print(f"step {i}\r", sep=' ', end='', flush=True)
    return [im]

writer = PillowWriter(fps=25)
anim = FuncAnimation(fig, animate, frames = 500, interval = 50)
anim.save('out.gif', writer=writer)
#plt.show()
