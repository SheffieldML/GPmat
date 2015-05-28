
Matlab Optimi Toolbox
=====================

This page describes the optimi toolbox for MATLAB available for [download here](http://www.cs.man.ac.uk/neill-bin/software/downloadForm.cgi?toolbox=optimi).
### Dependencies

The optimi toolbox depends on

- [ndlutil](https://github.com/SheffieldML/ndlutil) various utility functions.

### Release Information

Current release is 0.132.

Included a function to return the default functions to be used for particular constraints and a function to return the default optimiser.

Added wrappers for using Carl Rasmussen's conjugate gradient code.

This toolbox allows computation of several different non-linearities and associated conversions and their gradients.

There are several kernels implemented, the ones that are being maintained for the latest release are:

`exp sigmoid negLogLogit`

`exp` constrains a parameter to be positive through exponentiatiation. `sigmoid` constrains a parameter to be between 0 and 1 through the sigmoid function and `negLogLogit` constrains a parameter to be positive through the function `log(1+exp(x))`. In much of the code based on this we prefer `negLogLogit` to `exp` as it is better behaved for large `x`.

Page last modified on Fri Jan 5 22:19:06 GMT 2007.
