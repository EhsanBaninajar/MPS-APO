# MPS-APO
MPS-APO is a rapid and automatic parameter optimizer for multiple-point geostatistics

# Overview
MPS automatic parameter optimizer (MPS-APO), a generic method based on stochastic optimization to rapidly approximate optimal parameters for any MPS method and different kind of settings. The MPS-APO formulates an objective function that quantifies spatial pattern reproduction for each set of parameters. The Simultaneous Perturbation Stochastic Approximation (SPSA) optimization method is used because of its computational efficiency, and also its ability to cope with the stochastic nature of the objective function. The optimization proceeds in 2 steps. The first step aims to optimize the parameters for the best quality regardless of computational cost. When no more improvement can be achieved, the second step minimizes the CPU cost without degrading the spatial structures reproduction attained at the first step. MPS-APO is performed on different pixel-based and patch-based MPS methods: SNESIM, FILTERSIM, Direct Sampling and Image Quilting. The code is available for 2D and 3D multivariate Training Images for each method.

# Installation
Although any implementation on MPS methods can be used in MPS-APO, in this study for SNESIM and FILTERSIM the SGeMS implementation (Remy et al. 2009), for Direct Sampling the Mariethoz et al. (2010) implementation and for Image Quilting the Mahmud et al. (2014) implementation of the methods are used. 

### SNESIM and FILTERSIM
In order to used SNESIM-APO or FILTERSIM-APO, [mGstat](https://sourceforge.net/projects/mgstat/files/mGstat/) and [SGEMS](http://sgems.sourceforge.net/?q=node/77) integration in matlab is needed. 

### Direct Sampling
In order to acquire the DS code used in DS-APO, please contact Gregoire Mariethoz. 

### Image Quilting
The Matlab code for MPS simulation by Image Quilting is available [here](https://github.com/GAIA-UNIL/Image-Quilting). 

#
For any additional information or a bug report please contact me at: e.baninajar@gmail.com
