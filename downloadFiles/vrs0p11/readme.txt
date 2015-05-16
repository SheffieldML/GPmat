ReadMe file for the OPTIMI toolbox version 0.11 Thursday, June 17, 2004 at 15:43:44
Written by Neil D. Lawrence.

This toolbox implements functions which allow non-linear transformations between parameters to be optimised. For example it allows variances to be optimised in log space using expTransform. sigmoidTransform allows parameters to be constrained to be between 0 and 1, negLogLogitTransform constrains a parameter to be positive.

cmpndTieParameters.m: Tie parameters together.
expTransform.m: Constrains a parameter to be positive through exponentiation.
negLogLogitTransform.m: Constrains a parameter to be positive.
optimiseParams.m: Optimise parameters.
sigmoidTransform.m: Constrains a parameter to be between 0 and 1.
