#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

@np.vectorize
def f(t, N=10):
    return 1/N * int(N*t)

t = np.linspace(0, 1, 100)
for N,c in zip([2, 5, 10, 1000], ['b', 'r', 'g', 'black']):
    x = f(t, N=N)
    plt.figure()
    plt.step(t, x, 'o', linestyle='none', c=c, label='i = {}'.format(N))
    plt.grid()
    
plt.legend()
plt.show()
