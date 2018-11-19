import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

# Data for plotting
t = np.arange(0.0, 1.01, 0.01)
x = t*2-1
c = np.abs(x)

fig, ax = plt.subplots(figsize=plt.figaspect(.5))
plt.autoscale(tight=True)
ax.plot(x, c)

ax.set(xlabel='x', ylabel='y')
ax.grid()
plt.axis('scaled')
fig.savefig("figures/plot_abs_function.svg")
plt.show()
