ifgt-gpr
========

Build by running make bin/ifgt 

There are no special dependencies or libraries

Depending on which test is being run from main(), there is some default output. 

The IFGT tests output cluster information to out/clusters, and the 2d gauss transform on source location to out/2d
These files can be plotted using the python scripts in the plot/ folder by being passed as the primary argument. 
Plotting is done through matplotlib (on most Linux distros it is availabe under the package python-matplotlib).

Example workflow:
```
make bin/ifgt
./bin/ifgt
python plot/ifgt_imshow.py out/2d
python plot/kcenter.py out/clusters
```
This should produce the following plots: out/2d.pdf and out/clusters.pdf 
