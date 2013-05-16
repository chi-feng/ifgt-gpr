import matplotlib
matplotlib.use('PDF')

from matplotlib import cm, rc, pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
import os
from sys import argv, exit
argc = len(argv)

rc('text', usetex=True)
#rc('font', family='sans')
#rc('font', serif='Computer Modern Bright')
rc('font', size='12')
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{cmbright}"]

data = np.genfromtxt(argv[1])
x, y, z, e = data[:,0], data[:,1], data[:,2], abs(data[:,3])

xi = np.linspace(0.0, 1.0, 101)
yi = np.linspace(0.0, 1.0, 101)
X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z, xi, yi)
E = griddata(x, y, e, xi, yi)

fig = plt.figure()
fig.set_size_inches(6,4)
frame = plt.gca()

# plt.contourf(X, Y, Z, 100, cmap=cm.jet)
# levels = np.linspace(3.4, 3.7, 10)
# plt.contour(X, Y, Z, 15, cmap=cm.jet, linewidths=2)
#plt.title(r'$U(\mathbf{d})$ with  $\eta\sim\mathcal{N}(0,0.1)$')

# plt.clabel(cs, levels, inline=1, fmt='%1.1f', fontsize=14)
plt.subplot(1,2,1)
plt.imshow(Z, interpolation='none', origin='lower', cmap=cm.jet, extent=(0,1,0,1))
# plt.colorbar(format='%1f')
plt.subplot(1,2,2)
plt.imshow(E, interpolation='none', origin='lower', cmap=cm.jet, extent=(0,1,0,1))
plt.colorbar(format='%1.1e')

# set_backgroundcolor(plt.gca(), (0,0,0))

if argc > 2:
  plt.title(argv[2])

plt.savefig(argv[1] + '.pdf');
# plt.savefig(argv[1] + '.png', background='Transparent', dpi=150);
