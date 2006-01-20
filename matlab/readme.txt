The FGPLVM toolbox is a new implementation of the GP-LVM that uses the Pseudo-Input method of Snelson and Ghahramani (NIPS 2005) for sparsification and efficiency improvements. 

Version 0.13
------------

This version includes a much cleaner way of incorporating different dynamics models. It is released in line two imminent reports on learning large scale Gaussian processes and learning with back constraints.

Version 0.11
------------

This version now includes the Snelson-Ghahramani approximation (called FITC by Quinonero-Candela and Rasmussen) and the partially independent training criterion (PITC). Additionally the approximations can be used in standard Gaussian process regression.

Version 0.1
-----------

In the first release, only the projected latent variables approximation of Seeger et al is implemented. The toolbox also implements back-constraints as proposed by Lawrence and Quinonero-Candela

The first release containing a couple of examples on the oil data (demOil1.m and demOil2.m) and dynamics (demStick1.m and demStick2.m). The toolbox can load in the C++ files with dynamics associated and (through the mocapResultsCppBvh in the MOCAP toolbox) can run motion capture files with dynamics.