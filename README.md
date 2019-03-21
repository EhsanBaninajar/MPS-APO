# MPS-APO
MPS-APO is a rapid and automatic parameter optimizer for multiple-point geostatistics

# Overview
MPS automatic parameter optimizer (MPS-APO), a generic method based on stochastic optimization to rapidly approximate optimal parameters for any MPS method and different kind of settings. The MPS-APO formulates an objective function that quantifies spatial pattern reproduction for each set of parameters. The Simultaneous Perturbation Stochastic Approximation (SPSA) optimization method is used because of its computational efficiency, and also its ability to cope with the stochastic nature of the objective function. The optimization proceeds in 2 steps. The first step aims to optimize the parameters for the best quality regardless of computational cost. When no more improvement can be achieved, the second step minimizes the CPU cost without degrading the spatial structures reproduction attained at the first step. MPS-APO is performed on different pixel-based and patch-based MPS methods: SNESIM, FILTERSIM, Direct Sampling and Image Quilting. The code is available for 2D and 3D multivariate Training Images for each method.

# Installation
Any available MPS code can be used for the simulation:

SNESIM and FILTERSIM
In order to used SNESIM-APO or FILTERSIM-APO, Mgstat (https://sourceforge.net/projects/mgstat/files/mGstat/) and SGEMS (http://sgems.sourceforge.net/?q=node/77) integration in matlab is needed. 

