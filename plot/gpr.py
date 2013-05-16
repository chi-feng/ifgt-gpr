from matplotlib import cm, rc, pyplot as plt
import numpy as np
from sys import argv, exit

rc('text', usetex=True)
rc('font', size='12')
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{cmbright}"]

f = open(argv[1])

line = f.readline();
tokens = line.split(' ')
N = int(tokens[0])

x = []
y = []
noisefree = []
xgrid = []
mean = []
var = []
cg = []

for i in xrange(N):
  line = f.readline()
  tokens = line.split(' ')
  x.append(float(tokens[0]))
  y.append(float(tokens[1]))

line = f.readline();
tokens = line.split(' ')
M = int(tokens[0])

for i in xrange(M):
  line = f.readline()
  tokens = line.split(' ')
  xgrid.append(float(tokens[0]))
  mean.append(float(tokens[1]))
  var.append(float(tokens[2]))
  noisefree.append(float(tokens[3]))
  cg.append(float(tokens[4]))

mean = np.array(mean)
var = np.array(var)
cl = 1.96 * np.sqrt(var)
noisefree = np.array(noisefree)

plt.plot(x, y, 'k.', ms=10)
plt.plot(xgrid, noisefree, 'b--')
plt.plot(xgrid,mean,'k-',lw=2)
plt.plot(xgrid,cg,'k-',lw=1)
plt.fill_between(xgrid, mean+cl, mean-cl, facecolor='red', alpha=0.15)
plt.plot(xgrid,mean+cl,'r-',lw=1)
plt.plot(xgrid,mean-cl,'r-',lw=1)
plt.xlim([0,1])
plt.ylim([-1.2,1.2])
plt.savefig(argv[1] + '.pdf');
