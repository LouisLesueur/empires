import matplotlib.pyplot as plt
import numpy as np
from grid.grid import Grid
from maps.domain import Domain, Pop, Res

from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image


Europe = Domain('maps/europe')

N = 10
PI = []

count = 0
while(count) < N:
    i = np.random.randint(10, Europe.shape[0]//2-10)
    j = np.random.randint(10, Europe.shape[1]//2-10)

    if Europe.I[i,j] != 0:
        count +=1
        PI.append(Pop([i, j], 200, np.random.randint(255, size=3), 0.001, 10, 1000, 0.02, 0.0002, Europe.shape, Europe.area, Europe.dx))

Miam =  Res(0.2, 5000, 0.01, Europe.shape, Europe.area)

RHO = np.array([Miam])
sim = Grid(np.array(PI), RHO, Europe, 20)

fig = plt.figure(figsize=(10,5))

im=plt.imshow(sim.get_img())

def animate(i):
    for _ in range(10):
        sim.update()
    im.set_array(sim.get_img())

    print(f"step {i}\r", sep=' ', end='', flush=True)
    return [im]

writer = PillowWriter(fps=15)
anim = FuncAnimation(fig, animate, frames = 500, interval = 50)
#anim.save('out.gif', writer=writer)
plt.show()
