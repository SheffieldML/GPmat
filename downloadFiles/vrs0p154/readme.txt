NDLUTIL software
Version 0.154		Friday 05 Jan 2007 at 22:39
Copyright (c) 2007 Neil D. Lawrence

This toolbox implements some generic functions used by other toolboxes. Originally it was spun out of the IVM 0.221 toolbox.

Version 0.154
-------------

Added code for checking Hessian matrices.

Version 0.153
-------------

Added comments to files and changed jitChol, pdinv and logdet so that
they can return any jitter added. If the return the jitter then they don't emit a warning.

Version 0.152
-------------

There was a sign error in lnDiffCumGaussian. This has been fixed. The sign error was replicated in the NOISE and PRIOR toolboxes (a case of two wrongs making a 'right'). This means that earlier versions of this toolbox are incompatible with versions of PRIOR before 0.131 and versions of NOISE before 0.131.

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

jitChol.m: Do a Cholesky decomposition with jitter.
lnCumGaussian.m: log cumulative distribution for the normalised Gaussian.
ngaussian.m: Compute a Gaussian with mean 0 and variance 1.
xlogy.m: z = x*log(y) returns zero if x=y=0
traceProduct.m: Returns the trace of the product of two matrices.
kldivGaussian.m: Give the KL divergence between two Gaussians.
deg2rad.m: Transform degrees to radians.
gradientCheck.m: Check gradients of objective function.
scatterPlot.m: 2-D scatter plot of labelled points.
gaussSamp.m: Sample from a Gaussian with a given covariance.
defaultOptions.m: The default options for optimisation.
gaussOverDiffCumGaussian.m: A Gaussian over difference of cumulative Gaussians.
hessianCheck.m: Check Hessian of objective function.
invSigmoid.m: The inverse of the sigmoid function.
getSymbols.m: Get a cell array of different plot symbols.
lnDiffCumGaussian.m: Log of the difference between two cumulative Gaussians.
stack.m: Return column stacked vector of given matrix.
gradLogCumGaussian.m: Gradient of the log of the cumulative Gaussian.
tableRead.m: Read in data which has column titles in the first line and separated values in each other line.
getline.m: Get a line from a file.
numsf2str.m: Convert number to a string with a number of significant digits.
invCumGaussian.m: Computes inverse of the cumulative Gaussian.
cumGaussian.m: Cumulative distribution for Gaussian.
preparePlot.m: Helper function for tidying up the plot before printing.
negLogLogit.m: Function which returns the negative log of the logistic function.
rocCurve.m: Draw ROC curve and return labels.
logdet.m: The log of the determinant when argument is positive definite.
stringSplit.m: Return separate parts of a string.
pdinv.m: Invert a positive definite matrix.
zeroAxes.m: A function to move the axes crossing point to the origin.
tokenise.m: Split a string into separate tokens.
sigmoid.m: The sigmoid function
stringSigFigs.m: Convert number to a string with a number of significant digits.
sparseDiag.m: Create a diagonal matrix that is sparse from a vector.
