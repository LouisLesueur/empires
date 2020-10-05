import numpy as np


def grad(D, dx, axis=0):
    D_pad = np.pad(D, (1,1), 'reflect')
    return (np.roll(D_pad,1,axis=axis) - D_pad)[1:-1,1:-1]/dx

def grad_grad(D, A, dx):
    A_pad = np.pad(A, (1,1), 'reflect')
    D_pad = np.pad(D, (1,1), 'reflect')

    D_x = (np.roll(D_pad,1,axis=0) - D_pad)[1:-1,1:-1]
    A_x = (np.roll(A_pad,1,axis=0) - A_pad)[1:-1,1:-1]

    D_y = (np.roll(D_pad,1,axis=1) - D_pad)[1:-1,1:-1]
    A_y = (np.roll(A_pad,1,axis=1) - A_pad)[1:-1,1:-1]

    return (D_x*A_x + D_y*A_y)/(dx*dx)

def lap(A, dx):
    A_pad = np.pad(A, (1,1), 'reflect')
    return ((np.roll(A_pad,-1,axis=1) - 2*A_pad + np.roll(A_pad,1,axis=1)) + (np.roll(A_pad,-1,axis=0) - 2*A_pad + np.roll(A_pad,1,axis=0)))[1:-1,1:-1]/(dx*dx)
