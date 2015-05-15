NDLUTIL software
Version 0.151		Friday 24 Feb 2006 at 19:19
Copyright (c) 2006 Neil D. Lawrence

This toolbox implements some generic functions used by other toolboxes. Originally it was spun out of the IVM 0.221 toolbox.


Version 0.151
-------------

Sped up the sparseDiag command using spdiag, this makes the FGPLVM code run much faster.

Version 0.15
------------

Fixed inconsistencies in logdet and pdinv and corrected the amount of jitter that is added at first try. Placed the jitter addition in a new file jitChol which is called by both functions. 

Added kldivGaussian for computing Gaussian Kullback-Leibler divergences.

Added deg2rad for converting degrees to radians.

Version 0.142
-------------

Added sparseDiag and moved data loading methods into the new DATASETS toolbox. 

Version 0.14
------------

Moved in lvmLoadData and mappingLoadData as generic methods of loading in data sets. 

Version 0.131 Release Notes
---------------------------

Added files getline.m and tokenise.m which are string utilities for reading and processing files.

Version 0.13 Release Notes
--------------------------

pdinv now adds jitter which is proportional to the mean of the diagonal elements of the matrix.




MATLAB Files
------------

Matlab files associated with the toolbox are:

cumGaussian.m: Cumulative distribution for Gaussian.
defaultOptions.m: The default options for optimisation.
deg2rad.m: Transform degrees to radians.
gaussOverDiffCumGaussian.m: A gaussian in x over the difference between two cumulative Gaussians. 
gaussSamp.m: Sample from a Gaussian with a given covariance.
getline.m: Get a line from a file, but ignore it if it starts with a comment (default #).
getSymbols.m: Get a cell structure of different plot symbols.
gradientCheck.m: Check gradients of objective function.
gradLogCumGaussian.m: Gradient of the log of the cumulative Gaussian.
invCumGaussian.m: Inverser of the cumulative Gaussian.
invSigmoid.m: The inverse of the sigmoid function.
jitChol.m: Do a Cholesky decomposition, if matrix isn't positive definite add jitter and do it again.
kldivGaussian.m: Give the KL divergence between two Gaussians.
lnCumGaussian.m: log cumulative distribution for Gaussian.
lnDiffCumGaussian.m: Computes the log of the difference between two cumulative Gaussians.
logdet.m: The log of the determinant when argument is positive definite.
ngaussian.m: Compute a Gaussian with mean 0 and variance 1.
numsf2str.m: Convert number to a string with a number of significant digits.
pdinv.m: Invert a positive definite matrix.
preparePlot.m: Helper function for tidying up the plot before printing.
rocCurve.m: Draw ROC curve and return labels.
scatterPlot.m: 2-D scatter plot of labelled points.
sigmoid.m: The sigmoid function
sparseDiag.m: Create a diagonal matrix that is sparse from a vector.
stack.m: Return column stacked vector of given matrix.
stringSigFigs.m: Convert number to a string with a number of significant digits.
stringSplit.m: Return separate parts of a string.
tableRead.m: Read in data which has column titles in the first line and separated values in each other line.
tokenise.m: Split a string into separate tokens.
traceProduct.m: Returns the trace of the product of two matrices.
zeroAxes.m: A function to move the axes crossing point to the origin.
