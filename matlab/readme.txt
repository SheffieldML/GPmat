This toolbox allows priors to be placed over parameters, at the moment this is used so that MAP solutions can be found for the parameters rather than type-II maximum likelihood. The priors were written for the Null Category Noise Model (see NCNM toolbox) so that an exponential prior could be placed over the process variances. The rest of its functionality has not been tested, so use with care.


Version 0.132
-------------

Minor changes for reading in prior distributions from C++ written files.

Version 0.131
-------------

There was a sign error in lnDiffCumGaussian, and a corresponding sign error in the normal uniform prior, this has now been fixed.
