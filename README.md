[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/dektoud/Dual-Space-Clusterer/master?filepath=synthetic_example_agglomclus.ipynb)

# Dual-Space-Clusterer
Pure(ish) python implementation of the dual-space clustering algorithm.

Fortran modules are built using `compile.py` found in [`pygeostat`](http://www.ccgalberta.com/pygeostat/welcome.html)

The included `.pyd` files are built for windows 10 and python 3.6.3 64 bit with the intel fortran compiler and VS2015. 

NOTE: the included `*.ipynb` (and associated functions) requires the latest [`pygeostat`](http://www.ccgalberta.com/pygeostat/welcome.html)

This work is associated with the publication "Towards justifying unsupervised stationary decisions for geostatistical modeling: Ensemble spatial and multivariate clustering with geomodeling specific clustering metrics’" by Ryan Martin and Jeff Boisvert of the CCG at the University of Alberta, Edmonton, Canada
