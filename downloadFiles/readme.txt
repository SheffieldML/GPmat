OPTIMI software
Version 0.132		Friday 05 Jan 2007 at 22:20
Copyright (c) 2007 Neil D. Lawrence

The optimi toolbox is a helper toolbox that provides non-linear transformations for changing the regime across which parameters are optimised.

For example it allows variances to be optimised in log space using expTransform. sigmoidTransform allows parameters to be constrained to be between 0 and 1, negLogLogitTransform constrains a parameter to be positive.

Version 0.132
-------------

Included a function to return the default functions to be used for particular constraints and a function to return the default optimiser.

Version 0.131
-------------

Added wrappers for using Carl Rasmussen's conjugate gradient code.

Version 0.13
------------

Added optOptions as the default source for optimisation parameters. It checks if foptions is available, rather than assuming that it is.

MATLAB Files
------------

Matlab files associated with the toolbox are:

expTransform.m: Constrains a parameter to be positive through exponentiation.
cmpndTieParameters.m: Tie parameters together.
optimiseParams.m: Optimise parameters.
cgcarl.m: Wrapper for Carl Rasmussen's conjugate gradient implemntation.
optimiDefaultConstraint.m: Returns function for parameter constraint.
optimiDefaultOptimiser.m: Returns the default optimiser to be used.
gradFuncWrapper.m: Wrapper function to enable use of Carl Rasmussen's minimze function.
optOptions.m: Give optimisation options for NETLAB.
sigmoidTransform.m: Constrains a parameter to be between 0 and 1.
negLogLogitTransform.m: Constrains a parameter to be positive.
optimiMinimize.m: Wrapper for Carl Rasmussen's minimize function.
