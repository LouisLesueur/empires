import matplotlib.pyplot as plt
import numpy as np
from grid import Civ, Parcel, Simulation
from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image


I = plt.imread('maps/europe.png')
I_r = plt.imread('maps/europe_r.png')
I_R0 = plt.imread('maps/europe_r.png')

I = 0.2989 * I[:,:,0] + 0.5870 * I[:,:,1] + 0.1140 * I[:,:,2]
I_r = 0.2989 * I_r[:,:,0] + 0.5870 * I_r[:,:,1] + 0.1140 * I_r[:,:,2]
I_R0 = (0.2989 * I_R0[:,:,0] + 0.5870 * I_R0[:,:,1] + 0.1140 * I_R0[:,:,2])*100
I_Rmax = 5*I_R0

Tmax = 1000

sim = Simulation(I, I_r, I_R0, I_Rmax, Civ(0.8, 0.5, 1), [52, 66])
P = np.zeros((I.shape[0],I.shape[1]))
R = np.zeros((I.shape[0],I.shape[1]))

for i in range(I.shape[0]):
    for j in range(I.shape[1]):
        P[i,j] = sim.parcels[i][j].P
        R[i,j] = sim.parcels[i][j].R
sim.update()

fig = plt.figure()
ax = plt.axes()
im=plt.imshow(P,interpolation='none')


def animate(i):
    for i in range(I.shape[0]):
        for j in range(I.shape[1]):
            P[i,j] = sim.parcels[i][j].P
            R[i,j] = sim.parcels[i][j].R
    sim.update()

    im.set_array(P)
    return [im]

writer = PillowWriter(fps=25)
anim = FuncAnimation(fig, animate, frames = Tmax, interval = 0.00001)
anim.save('out.gif', writer=writer)
#plt.show()
