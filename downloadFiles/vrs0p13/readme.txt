ReadMe file for the NDLUTIL toolbox version 0.13 Monday, Jun 6, 2005 at 14:34:04
Written by Neil D. Lawrence.

License Info
------------

This software is free for academic use. Please contact Neil Lawrence if you are interested in using the software for commercial purposes.

This software must not be distributed or modified without prior permission of the author.


New Features:

pdinv now adds jitter which is proportional to the mean of the diagonal elements of the matrix.

This toolbox implements some generic functions used by other toolboxes. Originally it was spun out of the IVM 0.221 toolbox.

File Listing
------------

cumGaussian.m: Cumulative distribution for Gaussian.
defaultOptions.m: The default options for optimisation.
gaussOverDiffCumGaussian.m: A gaussian in x over the difference between two cumulative Gaussians. 
gaussSamp.m: Sample from a Gaussian with a given covariance.
getSymbols.m: Get a cell structure of different plot symbols.
gradLogCumGaussian.m: Gradient of the log of the cumulative Gaussian.
invCumGaussian.m: Inverser of the cumulative Gaussian.
invSigmoid.m: The inverse of the sigmoid function.
lnCumGaussian.m: log cumulative distribution for Gaussian.
lnDiffCumGaussian.m: Computes the log of the difference between two cumulative Gaussians.
logdet.m: The log of the determinant when argument is positive definite.
ndlutilVers.m: Brings dependent toolboxes into the path.
ngaussian.m: Compute a Gaussian with mean 0 and variance 1.
numsf2str.m: Convert number to a string with a number of significant digits.
pdinv.m: Computes the inverse of a positive definite matrix
preparePlot.m: Helper function for tidying up the plot before printing.
rocCurve.m: Draw ROC curve and return labels.
sigmoid.m: The sigmoid function
stringSigFigs.m: Convert number to a string with a number of significant digits.
stringSplit.m: Return separate parts of a string.
tableRead.m: Read in data which has column titles in the first line and separated values in each other line.
traceProduct.m: Returns the trace of the product of two matrices.
zeroAxes.m: A function to move the axes crossing point to the origin.
