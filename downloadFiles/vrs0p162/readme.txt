NDLUTIL software
Version 0.162		Friday 22 Jul 2011 at 16:00

This toolbox implements some generic functions used by other toolboxes. Originally it was spun out of the IVM 0.221 toolbox.

Version 0.162
-------------

Added files for printing to latex plots from Octave.

Version 0.161
-------------

Antti updated lnDiffErfs. Added readBinaryDoubles and printPlot.

Version 0.16
------------

Release for ICML tutorial. Added cumGamma and gammaPdf. Also added isoctave and several matlab helper files for reading from file IDs.

Version 0.159
-------------

Minor release for dimensional reduction demos. Added centeringMatrix.m.

Version 0.158
-------------

Antti Honkela provided the files lnDiffErfs and gradLnDiffErfs to assist in computing the DISIM kernel from the kernel toolbox stably.

Version 0.157
-------------

Added treeFindLeaves.

Version 0.156
-------------

Moved treeFindParents, treeFindChildren and treeSwapNode into this toolbox from MOCAP toolbox.

Version 0.155
-------------

Moved lnCumGaussSum from NCNM toolbox to this toolbox as part of merge of NCNM toolbox into NOISE and IVM toolboxes.


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

gammaPdf.m: PDF for the Gamma distribution.
preparePlot.m: Helper function for tidying up the plot before printing.
treeSwapNode.m: Swap two nodes in the tree structure array.
treeFindParents.m: Given a tree that lists only children, add parents.
cell2num.m: Converts a cell array of strings to numbers.
writeStringToFID.m: Writes a string to an FID.
readBoolFromFID.m: Read a boolean from an FID.
deg2rad.m: Transform degrees to radians.
invCumGaussian.m: Computes inverse of the cumulative Gaussian.
ngaussian.m: Compute a Gaussian with mean 0 and variance 1.
centeringMatrix.m: returns the centering matrix for the given dimensionality.
readStringFromFID.m: Read an boolean from an FID.
cumGaussian.m: Cumulative distribution for Gaussian.
xlogy.m: z = x*log(y) returns zero if x=y=0
safeSave.m: Safe save
sigmoid.m: The sigmoid function
gradLogCumGaussian.m: Gradient of the log of the cumulative Gaussian.
getSymbols.m: Get a cell array of different plot symbols.
cumGamma.m: Cumulative distribution for gamma.
jitChol.m: Do a Cholesky decomposition with jitter.
treeFindLeaves.m: Return indices of all leaf nodes in a tree structure.
lnCumGaussian.m: log cumulative distribution for the normalised Gaussian.
writeIntToFID.m: Writes an integer to an FID.
tokenise.m: Split a string into separate tokens.
treeGetWidths.m: give width of each level of tree.
tableRead.m: Read in data which has column titles in the first line and separated values in each other line.
readBinaryDoubles.m: Read information from a binary file in as doubles.
isoctave.m: Returns true if the software running is Octave.
writeBoolToFID.m: Writes a boolean to an FID.
scatterPlot.m: 2-D scatter plot of labelled points.
rocCurve.m: Draw ROC curve and return labels.
readDoubleFromFID.m: Read a double from an FID.
writeDoubleToFID.m: Writes a double to an FID.
traceProduct.m: Returns the trace of the product of two matrices.
printLatexText.m: Print a text string to file for latex input.
stringSplit.m: Return separate parts of a string.
xgamrnd.m: Draw a sample from the gamma distribution.
stack.m: Return column stacked vector of given matrix.
gradientCheck.m: Check gradients of objective function.
sparseDiag.m: Create a diagonal matrix that is sparse from a vector.
readIntFromFID.m: Read an integer from an FID.
pdinv.m: Invert a positive definite matrix.
escapeText.m: Add back slashes to escape existing backslashes in a text.
treeFindChildren.m: Given a tree that lists only parents, add children.
ndlutilToolboxes.m: Loads in toolboxes for NDLUTIL.
kldivGaussian.m: Give the KL divergence between two Gaussians.
treeFindRoots.m: Return indices of all root nodes in a tree structure.
printLatexPlot.m: Print a plot to LaTeX.
readVersionFromFID.m: Read version number from an FID.
stringSigFigs.m: Convert number to a string with a number of significant digits.
writeVersionToFID.m: Writes a version to an FID.
negLogLogit.m: Function which returns the negative log of the logistic function.
lnCumGaussSum.m: The log of the weighted sum of two cumulative Gaussians.
printPlot.m: Print a plot to eps and png files.
dist1.m:	Calculates absolute distance (i.e. L1 norm) between two sets of
invSigmoid.m: The inverse of the sigmoid function.
defaultOptions.m: The default options for optimisation.
gaussOverDiffCumGaussian.m: A Gaussian over difference of cumulative Gaussians.
numsf2str.m: Convert number to a string with a number of significant digits.
hessianCheck.m: Check Hessian of objective function.
gradLnDiffErfs.m: Compute the gradient of the log difference of two erfs.
getline.m: Get a line from a file.
logdet.m: The log of the determinant when argument is positive definite.
lnDiffCumGaussian.m: Log of the difference between two cumulative Gaussians.
zeroAxes.m: A function to move the axes crossing point to the origin.
gaussSamp.m: Sample from a Gaussian with a given covariance.
lnDiffErfs.m: Helper function for computing the log of difference
printLatexOptions.m: Options for printing a plot to LaTeX.
