import matplotlib.pyplot as plt
import numpy as np
from grid.grid import Grid
from maps.domain import Domain, Pop, Res, State

from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image


Europe = Domain('maps/europe')
Europe.resize(50)

lenN = 5
R =  Res(0.2,50)

count = 0
start = []
Nstart = []

while(count) < lenN:
    i = np.random.randint(10, Europe.shape[0]-10)
    j = np.random.randint(10, Europe.shape[1]-10)

    if Europe.I[i,j] != 0:
        count +=1
        start.append([i,j])
        Nstart.append(50)

N = Pop(start, #array with starting location coordonates
        Nstart, #start populations (pop)
        10, #maximum consumption rate per person (res^-1.pop^-1.year^-1)
        2, #R needed to reach half of maximum consumption rate (res)
        0.2, #pop natural growth rate (year^-1)
        2, #inflexion of KN
        10, #barbarian population level (pop)
        10) #migration inflexion factor (pop^-1)


sim = Grid(N, R, Europe, 0.01)
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
plt.colorbar()
plt.show()
