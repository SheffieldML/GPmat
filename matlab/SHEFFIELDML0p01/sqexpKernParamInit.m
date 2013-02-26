function kern = sqexpKernParamInit(kern)

% SQEXPKERNPARAMINIT SQEXP kernel parameter initialisation.
%
%	Description:
%	This kernel is a 'pre-packaged' compound kernel of the form
%	{'rbf', 'lin', 'bias', 'white'}. Using this kernel removes
%	the overhead of mutliple calls through the 'cmpnd' kernel.
%	
%	
%
%	KERN = SQEXPKERNPARAMINIT(KERN) initialises the pre-built compound
%	squared exponential kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	% SEEALSO SQEXPKERNPARAMINIT, KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2004 Neil D. Lawrence



kern.inverseWidth = 1;
kern.rbfVariance = 1;
kern.whiteVariance = 1; 
kern.biasVariance = 1;
kern.nParams = 4;

kern.transforms(1).index = [1 2 3 4];
kern.transforms(1).type = optimiDefaultConstraint('positive');

kern.isStationary = false;