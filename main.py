import matplotlib.pyplot as plt
import numpy as np
from grid.grid import Grid
from maps.domain import Domain, Pop, Res, State

from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image


Europe = Domain('maps/europe')
Europe.resize(15)

lenN = 5
R =  Res(0.3,600)

count = 0
start = []
Nstart = []

while(count) < lenN:
    i = np.random.randint(10, Europe.shape[0]//2-10)
    j = np.random.randint(10, Europe.shape[1]//2-10)

    if Europe.I[i,j] != 0:
        count +=1
        start.append([i,j])
        Nstart.append(500)

N = Pop(start, Nstart, 10, 4, 0.5, 1, 10, 10)

sim = Grid(N, R, Europe, 1)

fig = plt.figure()

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
