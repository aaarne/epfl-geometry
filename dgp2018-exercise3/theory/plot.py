#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

def f(t, N=10):
    return 1/N * int(N*t), t

t = np.linspace(0, 1, 100)
for N,c in zip([2, 5, 10], ['b', 'r', 'g']):
    for ti in t:
         plt.scatter(*f(ti, N=1000), c=c)

plt.plot(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
    
plt.grid()
plt.show()
