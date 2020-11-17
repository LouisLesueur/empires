import numpy as np


def lap(A, dx):
    A_pad = np.pad(A, (1,1), 'reflect')
    D_2_x = np.roll(A_pad,-1,axis=1) - 2*A_pad + np.roll(A_pad,1,axis=1)
    D_2_y = np.roll(A_pad,-1,axis=0) - 2*A_pad + np.roll(A_pad,1,axis=0)
    return (D_2_x+D_2_y)[1:-1,1:-1]/(dx*dx)

def bound_lap(A):
    A_pad = np.pad(A, (1,1), 'reflect')
    left = np.roll(A_pad,-1,axis=1)
    right = np.roll(A_pad,1,axis=1)
    down = np.roll(A_pad,-1,axis=0)
    up = np.roll(A_pad,1,axis=0)

    out = A_pad
    conflicts = np.zeros((*A_pad.shape,2))-1

    for v in [left,right,down,up,A_pad]:

        out[np.logical_and(out==-1,v==-1)] = -1
        out[np.logical_and(out==-1,v>-1)] = v[np.logical_and(out==-1,v>-1)]
        out[np.logical_and(np.logical_and(out>-1, v>-1), out != v)] = -2

    return out[1:-1,1:-1], conflicts[1:-1,1:-1]

def compute_bound(mask):
    bound_mask = np.zeros_like(mask)
    bound_mask[np.where(lap(mask,1)>0)] = 1
    bound_idx = np.array(np.where(bound_mask==1)).T
    bound_len = len(bound_idx)

    return bound_mask, bound_idx, bound_len
