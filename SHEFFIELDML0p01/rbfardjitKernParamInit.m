function kern = rbfardjitKernParamInit(kern)

% RBFARDJITKERNPARAMINIT RBFARD2 kernel parameter initialisation.
%
%	Description:
%	The automatic relevance determination version of the radial basis
%	function kernel (RBFARD2) is a very smooth non-linear kernel and is a
%	popular choice for generic use.
%	
%	k(x_i, x_j) = sigma2 * [exp(-1/2 *(x_i - x_j)'*A*(x_i - x_j)) + jit*delta_ij]
%	
%	The parameters are sigma2, the process variance (kern.variance), the
%	diagonal matrix of input scales (kern.inputScales, constrained to be
%	positive).
%	
%	
%
%	KERN = RBFARDJITKERNPARAMINIT(KERN) initialises the automatic
%	relevance determination radial basis function kernel structure with
%	some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%	
%
%	See also
%	RBFKERNPARAMINIT, KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias



% This parameter is restricted positive.
kern.variance = 1;
% this is a fixed parameter
kern.jitter = 1e-5;
kern.inputScales = 0.999*ones(1, kern.inputDimension);
kern.nParams = 1 + kern.inputDimension;

kern.transforms(1).index = [1:kern.nParams];
kern.transforms(1).type = optimiDefaultConstraint('positive');

kern.isStationary = true;