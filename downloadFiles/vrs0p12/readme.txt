ReadMe file for the OPTIMI toolbox version 0.12 Monday, July 12, 2004 at 19:25:18
Written by Neil D. Lawrence.

License Info
------------

This software is free for academic use. Please contact Neil Lawrence if you are interested in using the software for commercial purposes.

This software must not be distributed or modified without prior permission of the author.


This toolbox implements functions which allow non-linear transformations between parameters to be optimised. For example it allows variances to be optimised in log space using expTransform. sigmoidTransform allows parameters to be constrained to be between 0 and 1, negLogLogitTransform constrains a parameter to be positive.
File Listing
------------

cmpndTieParameters.m: Tie parameters together.
expTransform.m: Constrains a parameter to be positive through exponentiation.
negLogLogitTransform.m: Constrains a parameter to be positive.
optimiseParams.m: Optimise parameters.
sigmoidTransform.m: Constrains a parameter to be between 0 and 1.
