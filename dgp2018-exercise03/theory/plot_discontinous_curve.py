#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

def f(t):
    return t**2 * (t-1) * (t+1), t*(t-1)*(t+1)

t = np.linspace(0, 1, 100)
x = np.zeros((100, 2))
for i in range(100):
    x[i,0], x[i,1] = f(t[i])

plt.plot(x[:,0], x[:,1])
plt.show()
