Version 0.227
-------------

Added diag covariance which uses the inputs as a diagonal covariance function (takes a one dimensional input only). Useful for having heteroschedastic noise. And index covariance functions which return a value if the two input indices are identical, and zero otherwise. Also includes Jaakko's NDDISIM and NDSIM covariance functions.

Version 0.226
-------------

Added velotrans covariance function which allows a moving field type covariance function with a constant velocity. Added rbfperiodic2 for periodic covariances with learnable period. Added various multioutput covariance functions for switched latent force models.

Version 0.225
-------------

Updates for latest version of MULTIGP toolbox associated with MULTIGP tech report.

Version 0.224
-------------

Added disimSample and simSample for sampling from these multioutput covariance functions. Michalis added kernel types rbfard2 and linard2 which use a slightly different formulation of the ARD parameters.

Version 0.223
-------------

Minor fix of "ard" kernel which somehow had a kernel computation bit placed in the parameter initialization --- cut and past mistake at some point.

Version 0.222
-------------

Removed division by kernel variance in kernels for computing the variance of the kernel. It causes numerical problems when the variance is small. Also changed mlp kernels so that the default variance distant from the origin is 1 instead of pi/2. 

Version 0.221
-------------

Minor changes for reading in kernels
written by C++ code into files.

Version 0.22
------------

Added Wiener kernel and various kernels for multi output kernels
including white noise being propagated through the first and second
order differential equation. 


Version 0.21
------------

Compatibility changes for NCCA and SGPLVM toolboxes.

Version 0.2
-----------

Added multiKernGradX, added new kernels for diffusion processes and
2nd order differential equations.

Version 0.172
-------------

Further minor updates to kern for working with the new gpsim code
(corrected handling of white kernel in multiKern).


Version 0.171
-------------

Minor update to kernCreate to handle block kernel with differing
number of inputs for each block.

Version 0.17
------------

Further improvements on the stability of the sim kernel. Addition of
the driven input single input motif kernel (Antti Honkela) and the
modification of the multiKern type to allow each block to have a
different number of time points (Pei Gao).

Version 0.168
-------------

Found a bug in tensor gradient which meant gradients weren't being
computed correctly with respect to X when more X and X2 are both
provided as input arguments and both have length larger than 1.

Antti Honkela improved the numerial stability of the sim kernel
through judicious use of erfc.

Version 0.167
-------------

Added 'translate' kernel which allows wrapping of other kernels with a
kernel that translates the input location. Useful for moving the
non-stationarity around the input space.

Version 0.166
-------------

Added periodic version of RBF kernel (see Rasmussen and Williams pg 92
or Mackay's introduction to GPs (1998)) and periodic version of
Gibbs's non-stationary kernel (see e.g. pg 93 of Rasmussen and
Williams).

Version 0.165
-------------

Added flag which indicates whether or not a kernel is stationary. This
can be used for speeding computations (stationary kernels have a
constant diagonal). Also replaced calls to constraining functions with
'optimiDefaultConstraint' calls which return the default constraint
(making it easier for the user to change).

Version 0.163
-------------

This release removes the stubs for several KernDiagGradX.m files,
which were confusing kernDiagGradX.m, which assumes they only exist if
the function is implemented. For the kernel types 'lin', 'poly',
'mlp', and their 'ard' counter-types, these files existed but weren't
yet implemented.


Version 0.162
-------------

Added the Gibbs's non-stationary kernel, the rational quadratic kernel
and the Matern kernel with nu = 3/2 and nu = 5/2.

Version 0.161
-------------

Updated with the Single Input Motif kernel, and improved the
documentation.

Version 0.16
------------

An intermediate release with some problems.

Version 0.151
-------------
 
Added kernDiagGradient command and versions of it for rbf, white,
bias, whitefixed and rbfard code. This command improves the speed of
the fitc approximation in the FGPLVM code a lot.

Version 0.15
------------

Added tensor kernels and white noise kernels which don't return a
parameter for optimisation: 'whitefixed'. White fixed kernels were
added by Nathaniel King.

Version 0.142 Release Notes
---------------------------

Minor changes to the kernDisplay command and each command for the
relevant sub kernels. New command allows spaces to be placed in front
of the display so that the kernel can be better formated when
displayed as part of a larger model.

Also added is the kernel 'file' which is a kernel for which values are
precomputed and stored in a file.

Version 0.141 Release Notes
---------------------------

Included the kernGetVariance function for obtaining the `signal'
associated with the kernel for use with the FGPLVM toolbox.

Version 0.14 Release Notes
--------------------------

Added computation of parameter gradients with respect to sub-matrices
of the kernel matrix to allow for optimisation of inducing points.

Version 0.131 Release Notes
---------------------------

Added polynomial and polynomial ARD kernels for completeness (their
use is not recommended in Gaussian processes). Added kernReadFromFID.m
for reading in a kernel from a C++ written file.

Version 0.13 Release Notes
--------------------------

Added kernSetWhite as a helper function to set the level of white
noise in the kernel.

General Overview
----------------

This toolbox implements the different kernels. At the time of writing
two toolboxes make use of KERN, IVM vs 0.31 and FBD vs 0.2.

Interaction with the toolbox is done through the interface files which
are prefixed by kern. The toolbox is designed to allow linear
combinations of kernels with a minimum of fuss (using the cmpnd
kernel).

The toolbox was spun out of the IVM toolbox, and most of the files are
based on files in IVM 0.221.

Kernel Types
------------

Several example kernels are given:

     'ard' For backward compatability with the ard kernel in IVM 0.1. It combines linear and rbf ard kernels.

     'sqexp' For backwards compatability this is equivalent to the 'rbf' kernel in IVM 0.1.

     'mlp', 'mlpard' The multi-layer perceptron kernel from Williams' Computing with infinite networks paper. An ARD version is also provided.

     'rbf', 'rbfard' The standard radial basis function kernel and an ARD version.

     'lin', 'linard' A linear kernel and an ARD version.

     'white' Is just a white noise kernel. It is not designed to be used alone, but as an element in the compound kernel

     'bias' is for adding a bias variance term to the kernel (a positive offset) on it's own it is not a valid kernel.

     'cmpnd' The compound kernel is for creating new kernels which are linear combinations of other kernels.

The perl script for generating code for new kernels is kernelGenerator.pl

It is run with two arguments, the first is the short name for the noise model, e.g. rbf, the second is the long name, e.g. radial\ basis\ function.
