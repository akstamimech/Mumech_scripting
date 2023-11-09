import subprocess, os
import numpy as np 
import matplotlib.pyplot as plt 
import csv
from pathlib import Path
from tqdm import tqdm
import pathlib
import warnings 
from mpl_toolkits.mplot3d import axes3d
# from mayavi import mlab
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
warnings.filterwarnings("ignore")


pointsx = [104,104,104,104,52,52,52,13,13,13,13,7,7,7]
pointsy = [0.0478,0.1912,0.5736,0.3824,0.0478,0.3824,0.5736,0.0956,0.1912,0.3824,0.5736,0.1912,0.3824,0.5736]
pointsz = [0.001281,0.0009856,0.0007029,0.0006898,0.001106,0.0007006,0.0006362,0.0005908,0.0006014,0.0006062,0.000598,0.000582,0.0005955,0.0005439]


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')



ax.set_xticks(pointsx)
ax.set_yticks(pointsy)
# ax.set_zticks(pointsz)

ax.plot_trisurf(pointsx, pointsy, pointsz, cmap=cm.jet, linewidth=0, alpha = 0.5)
ax.plot([104, 52], [0.3824, 0.3824], [0.0006898,0.0007006], color = "red")
# fig.colorbar(surf)
ax.scatter(pointsx, pointsy, pointsz)
ax.set_xlabel("FG spacing")
ax.set_ylabel("c/h ratio")
ax.set_zlabel("Diffusion coefficient")


# ax.xaxis.set_major_locator(MaxNLocator(5))
# ax.yaxis.set_major_locator(MaxNLocator(6))
# ax.zaxis.set_major_locator(MaxNLocator(5))



fig.tight_layout()

# plt.show() # or:
# fig.savefig('3D.png')

# ax = plt.figure().add_subplot(projection='3d')

# Plot the 3D surface
# ax.plot(pointsx, pointsy, pointsz, lw=0.5,
                # alpha=0.3)

# Plot projections of the plots for each dimension.  By choosing offsets
# that match the appropriate axes limits, the projected plots will sit on
# the 'walls' of the graph.
# ax.plot(pointsx, pointsy, zdir='z')
# ax.plot([0.0478, 0.1912,0.3824,0.5736],[0.001281,0.0009856,0.0006898,0.0007029] , zdir='x', label = "fg104")
# ax.plot([0.0478,0.3824,0.5735], [0.001106,0.0007006,0.0006362], zdir = 'x', label = "fg52")
# ax.plot([0.0956,0.1912,0.3824,0.5736], [0.0005908,0.0006014,0.0006062,0.000598], zdir = "x", label = "fg13")
# ax.plot([0.1912,0.3824,0.5736], [0.000582,0.0005955,0.0005439], zdir = "x", label = "fg7")
# ax.scatter(pointsx, pointsz, zdir='y')

ax.axes.set_zlim3d(bottom=0.0005)
# ax.set(xlim=(-40, 40), ylim=(-40, 40), zlim=(-100, 100),
#        xlabel='X', ylabel='Y', zlabel='Z')

plt.show()