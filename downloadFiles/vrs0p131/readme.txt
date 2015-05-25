OPTIMI software
Version 0.131		Saturday 14 Oct 2006 at 23:01
Copyright (c) 2006 Neil D. Lawrence



The optimi toolbox is a helper toolbox that provides non-linear transformations for changing the regime across which parameters are optimised.

For example it allows variances to be optimised in log space using expTransform. sigmoidTransform allows parameters to be constrained to be between 0 and 1, negLogLogitTransform constrains a parameter to be positive.

Version 0.131
-------------

Added wrappers for using Carl Rasmussen's conjugate gradient code.

Version 0.13
------------

Added optOptions as the default source for optimisation parameters. It checks if foptions is available, rather than assuming that it is.

MATLAB Files
------------

Matlab files associated with the toolbox are:

cgcarl.m: Wrapper for Carl Rasmussen's conjugate gradient implemntation.
cmpndTieParameters.m: Tie parameters together.
expTransform.m: Constrains a parameter to be positive through exponentiation.
gradFuncWrapper.m: Wrapper function to enable use of Carl Rasmussen's minimze function.
negLogLogitTransform.m: Constrains a parameter to be positive.
optimiseParams.m: Optimise parameters.
optOptions.m: Give optimisation options for NETLAB.
sigmoidTransform.m: Constrains a parameter to be between 0 and 1.
