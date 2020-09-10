import matplotlib.pyplot as plt
import numpy as np
from simulation.civ import Civ
from simulation.parcel import Parcel
from simulation.grid import Grid
from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image


I = plt.imread('maps/europe.png')
I_r = plt.imread('maps/europe.png')
I_R0 = plt.imread('maps/europe.png')

I = 0.2989 * I[:,:,0] + 0.5870 * I[:,:,1] + 0.1140 * I[:,:,2]
I_r = (0.2989 * I_r[:,:,0] + 0.5870 * I_r[:,:,1] + 0.1140 * I_r[:,:,2])
I_R0 = (0.2989 * I_R0[:,:,0] + 0.5870 * I_R0[:,:,1] + 0.1140 * I_R0[:,:,2])*12

I_r *= 0.08

I_Rmax = 5*I_R0

Tmax = 500

Rome = [50, 62]
Romans = Civ(0.01, 0.03, 0.02, 5, Rome, np.array([1,0,0]))

Barbaria = [20, 120]
Barbarians = Civ(0.02, 0.03, 0.02, 1, Barbaria, np.array([1,0,1]))


sim = Grid(I, I_r, I_R0, I_Rmax, [Romans, Barbarians] )

#P = np.zeros((I.shape[0],I.shape[1]))
#R = np.zeros((I.shape[0],I.shape[1]))
Col = np.zeros((I.shape[0],I.shape[1],3))

for i in range(I.shape[0]):
    for j in range(I.shape[1]):
        #P[i,j] = sim.parcels[i][j].P
        #R[i,j] = sim.parcels[i][j].R
        Col[i,j] = sim.parcels[i][j].civ.color
sim.update()

fig = plt.figure()
ax = plt.axes()
im=plt.imshow(Col,interpolation='none')


def animate(i):
    for i in range(I.shape[0]):
        for j in range(I.shape[1]):
            Col[i,j] = sim.parcels[i][j].civ.color
    sim.update()

    im.set_array(Col)
    return [im]

writer = PillowWriter(fps=25)
anim = FuncAnimation(fig, animate, frames = Tmax, interval = 0.00001)
#anim.save('out.gif', writer=writer)
plt.show()
