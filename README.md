Matlab Kernel Toolbox
=====================


Release Information
-------------------

**Current release is 0.227**.

As well as downloading the KERN software you need to obtain the toolboxes specified below. 

- [optimi](https://github.com/SheffieldML/optimi) optimization constriant mappings.
- [ndlutil](https://github.com/SheffieldML/ndlutil) various utility functions.

- [netlab](https://github.com/sods/netlab) Ian Nabney's Netlab toolbox.

- [erfcxz](http://www.mathworks.com/matlabcentral/fileexchange/12091-complex-scaled-complementary-error-function) Thomas Winiecki's erfcxz function. You need to rename W() to erfcxz to ensure it runs. It is available from the MATLAB central file exchange.

- [erfz](http://www.mathworks.com/matlabcentral/fileexchange/3574-erfz) Paul Godfrey's erfz Function. It is available from the MATLAB central file exchange.

Added diag covariance which uses the inputs as a diagonal covariance function (takes a one dimensional input only). Useful for having heteroschedastic noise. And index covariance functions which return a value if the two input indices are identical, and zero otherwise. Also includes Jaakko's NDDISIM and NDSIM covariance functions.

### Version 0.226

Added velotrans covariance function which allows a moving field type covariance function with a constant velocity. Added rbfperiodic2 for periodic covariances with learnable period. Added various multioutput covariance functions for switched latent force models.

### Version 0.225

Updates from Mauricio for the latest release of the MULTIGP toolbox.

### Version 0.224

Added disimSample and simSample for sampling from these multioutput covariance functions. Michalis added kernel types rbfard2 and linard2 which use a slightly different formulation of the ARD parameters.

### Version 0.223

Minor fix of "ard" kernel which somehow had a kernel computation bit placed in the parameter initialization --- cut and past mistake at some point.

### Version 0.222

Removed division by kernel variance in kernels for computing the variance of the kernel. It causes numerical problems when the variance is small. Also changed mlp kernels so that the default variance distant from the origin is 1 instead of pi/2.

### Version 0.221

Fixed code for reading in kernels from C++ files.

### Version 0.22

Added Wiener kernel and various kernels for multi output kernels including white noise being propagated through the first and second order differential equation.

### Version 0.21

Compatibility changes for NCCA and SGPLVM toolboxes.

### Version 0.2

Further minor updates to kern for working with the new gpsim code (corrected handling of white kernel in multiKern).

### Version 0.171

Minor changes to kernCreate for multiKern structures where there are different numbers of points in each block.

### Version 0.17

Further improvements on the stability of the sim kernel. Addition of the driven input single input motif kernel (Antti Honkela) and the modification of the multiKern type to allow each block to have a different number of time points (Pei Gao).

### Version 0.168

Found a bug in tensor gradient which meant gradients weren't being computed correctly with respect to X when more X and X2 are both provided as input arguments and both have length larger than 1.

Antti Honkela improved the numerial stability of the sim kernel through judicious use of erfc.

### Version 0.167

Added 'translate' kernel which allows wrapping of other kernels with a kernel that translates the input location. Useful for moving the non-stationarity around the input space.

### Version 0.166

Added periodic version of RBF kernel (see Rasmussen and Williams pg 92 or Mackay's introduction to GPs (1998)) and periodic version of Gibbs's non-statinary kernel (see e.g. pg 93 of Rasmussen and Williams).

### Version 0.165

Added flag which indicates whether or not a kernel is stationary. This can be used for speeding computations (stationary kernels have a constant diagonal). Also replaced calls to constraining functions with 'optimiDefaultConstraint' calls which return the default constraint (making it easier for the user to change).

### Version 0.163

This release removes the stubs for several KernDiagGradX.m files, which were confusing kernDiagGradX.m, which assumes they only exist if the function is implemented. For the kernel types 'lin', 'poly', 'mlp', and their 'ard' counter-types, these files existed but weren't yet implemented.

### Version 0.162

Added the Gibbs's non-stationary kernel, the rational quadratic kernel and the Matern kernel with nu = 3/2 and nu = 5/2.

### Verison 0.161

Introduced the single input motif kernel for the GPSIM toolbox. Also there is much more documentation, and a new file kernelGenerator.py for creating the basic files for your own kernels.

Examples
--------

This toolbox allows computation of several different kernels and their gradients. You can add kernels to the toolbox by creating versions of the relevant files. Once added, they can be tested using the `kernTest`. For example you can test the RBF kernel by writing

```matlab
>> kernTest('rbf')
```

There are several kernels implemented, the ones that are being maintained for the latest release are:

`  gibbs  gibbsperiodic  lin  linard  rbf  rbfard  rbfperiodic  matern32  matern52  ratquad  mlp  mlpard  poly  polyard  sim  lfm  disim  white  whitefixed  bias  cmpnd  wiener  gg  ou  lfmwhite  simwhite  ggwhite  gaussianwhite  gaussian  tensor (tensor kernels).  file (a kernel written in a file). `

A new kernel can be created using the compound kernel, `cmpnd` or the tensor kernel, `tensor`. The compound kernel is made up of sums of individual kernels. Most often you will need to create a kernel containing a base kernel (e.g. `rbf`) and adding a white noise kernel, `white` and perhaps a constant offset through the bias kernel `bias`. On initialisation most kernel parameters are set to 1. Exceptions are ARD scale parameters, the variance of the white and bias kernels and the weight and bias variances of the `mlp` kernels.

LFM Kernel
----------

### FORTRAN Compiler for LFM Kernel

To install the compiler

Go to <http://www.g95.org/downloads.shtml> and download the binary version suitable for your computer and operating system or compile the source code.

In MATLAB, write

```matlab
>>  mex -setup
```

Then choose option 1 and make a copy `mexopts.sh` in your local directory. Change the name to `g95opts.sh` or whatever you prefer.

Modify myopts.sh following the instructions in <http://www.g95.org/howto.shtml#matlab>

When compiling in MATLAB, use the -f command to use your local `g95opts.sh` file, for example:

```
>> mex -f myopts.sh lfmComputeUpsilonMatrix.f
```

Page updated on Tue Aug 9 20:39:05 2011