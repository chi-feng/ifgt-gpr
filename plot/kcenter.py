import matplotlib
from matplotlib import cm, rc, pyplot as plt
import numpy as np
from numpy.random import rand
import os
from sys import argv, exit
argc = len(argv)

rc('text', usetex=True)
rc('font', size='12')
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{cmbright}"]


def hsvToRGB(h, s, v):
    hi = np.floor(h / 60.0) % 6
    f =  (h / 60.0) - np.floor(h / 60.0)
    p = v * (1.0 - s)
    q = v * (1.0 - (f*s))
    t = v * (1.0 - ((1.0 - f) * s))
    return {
        0: (v, t, p),
        1: (q, v, p),
        2: (p, v, t),
        3: (p, q, v),
        4: (t, p, v),
        5: (v, p, q),
    }[hi]

fig = plt.figure()
fig.set_size_inches(5,4)
frame = plt.gca()

f = open(argv[1])
line = f.readline();
tokens = line.split(' ')
K = int(tokens[0])
N = int(tokens[1])

centers = []

groups = []

for i in xrange(K):
  line = f.readline()
  tokens = line.split(' ')
  centers.append(int(tokens[1]))
  groups.append([])

points = []

for i in xrange(N):
  line = f.readline()
  tokens = line.split(' ')
  points.append([int(tokens[1]), float(tokens[2]), float(tokens[3])])
  groups[int(tokens[1])].append([float(tokens[2]), float(tokens[3])])
  
for group in groups:
  x, y = zip(*group)
  color = hsvToRGB(rand() * 255, 1.0, rand() * 0.5 + 0.5)
  plt.plot(x,y,'.',color=color, mew=0, ms=2)
  
print len(centers)

for center in centers:
  plt.plot(points[center][1],points[center][2],'x',color='k',mew=0.5,ms=3)
  
plt.savefig(argv[1] + '.png', dpi=150);
