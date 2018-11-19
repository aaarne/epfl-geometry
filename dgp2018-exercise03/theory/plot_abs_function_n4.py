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

lines = [ mlines.Line2D([-1,-0.5],[1,0.5], color='r', linestyle='-', marker='x'),\
			mlines.Line2D([-0.5,0.5],[0.5,0.5], color='r', linestyle='-'),\
			mlines.Line2D([0.5,1],[0.5,1], color='r', linestyle='-', marker='x')]

for l in lines:
	ax.add_line(l)

ax.set(xlabel='x', ylabel='y')
ax.grid()
plt.axis('scaled')

fig.savefig("figures/plot_abs_function_n4.svg")
plt.show()
