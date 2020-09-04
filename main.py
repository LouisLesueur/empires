import matplotlib.pyplot as plt
import numpy as np
from grid import Civ, Parcel, Grid
import matplotlib.animation as animation

test_civ = Civ(0.5, 0.5, 1)

test_grid = Grid(30, 30, [15,15])

Tmax = 1000


P = np.zeros((30,30,Tmax))
R = np.zeros((30,30,Tmax))

for t in range(Tmax):
    for i in range(30):
        for j in range(30):
            P[i,j,t] = test_grid.parcels[i][j].P
            R[i,j,t] = test_grid.parcels[i][j].R



    test_grid.update()


plt.imshow(R[:,:,0])
plt.colorbar()
plt.plot()

fig = plt.figure()
ax = plt.axes(xlim=(0, 30), ylim=(0, 30))
im=plt.imshow(P[:,:,0],interpolation='none')


def animate(i):

    im.set_array(P[:,:,i])
    return [im]


anim = animation.FuncAnimation(
                               fig,
                               animate,
                               frames = Tmax,
                               interval = 1000 / 24, # in ms
                               )

plt.show()
