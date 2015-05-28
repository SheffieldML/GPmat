Changes in vs 0.41
------------------

Fixed bug where, for multiple outputs, the IVM 'random' point inclusion wasn't computing the entropy change in ivmSelectPoint.


Changes in vs 0.4
-----------------

The toolbox has been brought into line with other toolboxes on the site. Now the ivm model is created using ivmCreate, and all the model options are set in ivmOptions. The files are now better descrived. The null category noise model examples are now integrated with the IVM (the main ncnm code is in the NOISE toolbox).

Changes in vs 0.33
------------------

The code now relies upon the datasets toolbox for loading of datasets used. The EP updates have now been fixed so that the site parameters can be refined.


Changes in vs 0.32
------------------

This code no longer works under MATLAB 6.1, it requires MATLAB 7.0 or higher to run.
 
This version of the software now relies on the following toolboxes:

KERN vs 0.14
------------

This toolbox implements the different kernels. IVM interacts with this toolbox through an interface which involves files starting with kern.

NOISE vs 0.12
-------------

This toolbox implements the different noise models. IVM interacts with this toolbox through an interface which involves files starting with noise.

NDLUTIL vs 0.12
---------------

This toolbox implements some generic functions which could be used beyond the ivm toolbox, for example sigmoid.m, cumGaussian.m

OPTIMI vs 0.12
--------------

This toolbox implements functions which allow non-linear transformations between parameters to be optimised. For example it allows variances to be optimised in log space.

PRIOR vs 0.12
-------------

This toolbox allows priors to be placed over parameters, at the moment this is used so that MAP solutions can be found for the parameters rather than type-II maximum likelihood. The priors were written for the Null Category Noise Model (see NCNM toolbox) so that an exponential prior could be placed over the process variances. The rest of its funcitonality has not been tested.

ROCHOL vs 0.12
--------------

This toolbox implements the rank one Cholesky updates. These are need when points are removed or EP updates are applied to selected points.

Changes in vs 0.31
------------------

The options are now held in a structure whose values are set in ivmOptions.m
There were some missing files in the last release, these have now been added.

The EP updates are currently unstable numerically and should be used with caution.

Demos
-----

Six toy demos are provided: demClassification1, demClassification2, demRegression1, demRegression2, demOrdered1 and demOrdered2. Each runs a different noise model with. They display their results as they run and therefore they don't use the ivmRun function which is the normal recommended way for running code.

Three large scale experiments are provided on the USPS data-set, demUsps1-3. The use three different types of kernel.

General comments
----------------

Since version 0.22 the code is far more modular, this was done in an effort to improve its readability and reduce the need for re-writes. However it may be slower than the previous version as a result. 

Yet to be implemented functionality still includes:

Multi-class noise models (which will probably be done mostly in the NOISE toolbox) and randomised greed selection.
