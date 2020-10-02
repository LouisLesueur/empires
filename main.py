import matplotlib.pyplot as plt
import numpy as np
from simulation.grid import Grid

from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image


I = plt.imread('maps/europe.png')
I_r = plt.imread('maps/europe_r.png')
I_R0 = plt.imread('maps/europe_r.png')

I = 0.2989 * I[:,:,0] + 0.5870 * I[:,:,1] + 0.1140 * I[:,:,2]
I_r = (0.2989 * I_r[:,:,0] + 0.5870 * I_r[:,:,1] + 0.1140 * I_r[:,:,2])*0.0005
I_R0 = (0.2989 * I_R0[:,:,0] + 0.5870 * I_R0[:,:,1] + 0.1140 * I_R0[:,:,2])*100+1

N1 = np.zeros_like(I)
N2 = np.zeros_like(I)
N3 = np.zeros_like(I)
N4 = np.zeros_like(I)

N1[50, 62] = 1  #Rome
N2[30, 83] = 1 #Barbarians
N3[34, 44] = 1  #Paris
N4[10, 134] = 1  #Russe


sim = Grid(np.array([I_R0]), #R
           np.array([N1,N2,N3,N4]),
           0.005*np.array([np.ones(N1.shape), np.ones(N2.shape), np.ones(N2.shape), np.ones(N2.shape)]), #gamma
           np.array([I_r]), #r
           np.array([I_R0]), #KR
           0.05*np.array([1, 1, 1, 1]), #DN_0
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

def colormap(pi, color, max):
    out = color*np.ones((pi.shape[0], pi.shape[1], 3))
    fact = pi/max
    #fact[np.where(pi>0.001)] = 1


    out[:,:,0] *= fact
    out[:,:,1] *= fact
    out[:,:,2] *= fact
    return out


I = plt.imread('maps/europe.png')

colors = np.array([[0,0,1],[0,1,0],[1,0,0], [1,0,1]])
out = np.zeros((sim.pi[0].shape[0], sim.pi[0].shape[1], 3))
for i in range(len(sim.pi)):
    out += colormap(sim.pi[i], colors[i], 1)

fig = plt.figure(figsize=(10,5))
fig.tight_layout(pad=10) # Or equivalently,  "plt.tight_layout()"

ax1 = fig.add_subplot(1, 2, 1)
im=plt.imshow((out*255).astype(np.uint8))

ax2 = fig.add_subplot(1, 2, 2)
im2, = plt.plot([],[],label="pop1",color = colors[0])
im21, = plt.plot([],[],label="pop2",color = colors[1])
im22, = plt.plot([],[],label="pop3",color = colors[2])
im23, = plt.plot([],[],label="pop4",color = colors[3])
im24, = plt.plot([],[],label="ressources")
plt.legend()

plt.ylim(0,1)
out2 = []
out21 = []
out22 = []
out23 = []
out24 = []


def animate(i):
    out = np.zeros((sim.pi[0].shape[0], sim.pi[0].shape[1], 3))
    for k in range(sim.pi.shape[0]):
        out += colormap(sim.pi[k], colors[k], np.max(sim.pi))
    for _ in range(10):
        sim.update()
    im.set_array((out*255).astype(np.uint8))
    out2.append(np.sum(sim.pi[0]))
    out21.append(np.sum(sim.pi[1]))
    out22.append(np.sum(sim.pi[2]))
    out23.append(np.sum(sim.pi[3]))
    out24.append(np.sum(sim.rho[0]))

    print(np.min(sim.a[1,0]), np.max(sim.a[1,0]), np.mean(sim.a[1,0]))

    im2.set_data(np.arange(i+2), out2/np.max(out2))
    im21.set_data(np.arange(i+2), out21/np.max(out21))
    im22.set_data(np.arange(i+2), out22/np.max(out22))
    im23.set_data(np.arange(i+2), out23/np.max(out23))
    im24.set_data(np.arange(i+2), out24/np.max(out24))

    ax2.set_xlim(0,i)
    print(f"step {i}\r", sep=' ', end='', flush=True)
    return [im]

writer = PillowWriter(fps=15)
anim = FuncAnimation(fig, animate, frames = 10000, interval = 50)
#anim.save('out.gif', writer=writer)
plt.show()
