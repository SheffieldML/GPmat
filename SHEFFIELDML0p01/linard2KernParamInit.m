function kern = linard2KernParamInit(kern)

% LINARD2KERNPARAMINIT LINARD2 kernel parameter initialisation.
%
%	Description:
%	The automatic relevance determination version of the linear
%	kernel (LINARD2) is the simple inner product kernel with feature
%	selection applied.
%	
%	k(x_i, x_j) = x_i'*A* x_j
%	
%	where A is a diagonal matrix of values constrained to positve. These
%	parameters are stored in the field 'inputScales'.
%	
%	
%
%	KERN = LINARD2KERNPARAMINIT(KERN) initialises the automatic
%	relevance determination linear kernel structure with some default
%	parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%	
%
%	See also
%	LINKERNPARAMINIT, RBFARDKERNPARAMINIT, KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias


% These parameters are restricted to positive
kern.inputScales = 0.999*ones(1, kern.inputDimension);
kern.nParams = kern.inputDimension;

kern.transforms(1).index = [1:kern.nParams];
kern.transforms(1).type = optimiDefaultConstraint('positive');

kern.isStationary = false;