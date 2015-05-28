Version 0.141 Release Notes
---------------------------

Added noiseReadFromFile and noiseReadFromFID for compatability with CPP releases.

Version 0.14 Release Notes
--------------------------

Much improved documentation and merging of NCNM noise model into this toolbox from the NCNM toolbox.

Version 0.131 Release Notes
---------------------------

There was a sign error in lnDiffCumGaussian and a corresponding sign error in orderedLogLikelihood. This has now been fixed.


Version 0.13 Release Notes
--------------------------

This code no longer runs under MATLAB 6.0. It now requires MATLAB 7.0 or higher to run. 

Version 0.121 Release Notes
---------------------------

Included a `scale' noise type for use with the GPLVM which allows for scaling the outputs on the GPLVM. This noise model is not designed for normal use with the IVM.

Added noiseReadFromFID.m for reading a noise model from a file written by the C++ code.

Version 0.12 Release Notes
--------------------------

This toolbox implements different noise models for the IVM toolbox from version 0.31.

Interaction with the toolbox is done through the interface files which are prefixed by `noise'. 

The toolbox was spun out of the IVM toolbox, and most of the files are based on files in IVM 0.221.

Noise Types
-----------

Three main noise models are also provided:

      'gaussian' is the standard Gaussian noise model.

      'probit' is the probit model for classification.

      'ordered' is an ordered categorical model based on the probit.

Also there is

      'cmpnd' is for associating different noise models to different processes when learning multiple processes together. The ability to learn multiple processes is mainly included so that the next release of GPLVM code can use this IVM code base, but it may also be useful for multi-class classification noise models.

      'mgaussian' which is designed to allow multiple processes to have individually different variances (mainly for the GPLVM).

The perl script for generating code for new noise models is noiseGenerator.pl

It is run with two arguments, the first is the short name for the noise model, e.g. gaussian, the second is the long name, e.g. Gaussian\ noise\ model.
