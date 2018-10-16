#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 0.5*np.pi, 100)
y = np.sin(x)*2*np.pi
plt.plot(x,y)
plt.grid()
plt.xlabel('parameter v')
plt.ylabel('area density')
plt.savefig('figures/density.svg')
