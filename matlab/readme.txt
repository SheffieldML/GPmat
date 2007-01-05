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