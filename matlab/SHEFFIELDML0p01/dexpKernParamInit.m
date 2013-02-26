function kern = dexpKernParamInit(kern)

% DEXPKERNPARAMINIT The double exponential kernel is usually called
%
%	Description:
%	simply the exponential kernel in the machine learning literature, and
%	sometimes also the Laplace kernel. This name is chosen because the
%	expression of the kernel,
%	
%	k(x_i, x_j) = 0.5 * sigma2 * theta * exp(-theta*abs(x_i - x_j)),
%	
%	is identical to that of the multivariate Laplace probability density
%	function (pdf), which is also called sometimes double exponential pdf
%	(at least in 1D).
%	
%	The parameters are sigma2, the process variance (kern.variance), and
%	theta, the inverse width of the kernel (kern.inverseWidth).
%	
%	This is a stationary kernel and its one-dimensional version is similar
%	to the stationary version of the OU kernel. However, unlike the OU kernel
%	there is no constraint in the range of the inputs.
%	
%	
%
%	KERN = DEXPKERNPARAMINIT(KERN) initialises the double exponential
%	kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%	
%
%	See also
%	OUKERNPARAMINIT, KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2009 David Luengo
%	Copyright (c) 2009 Neil D. Lawrence

  


% Parameters of the cross-covariance function
kern.nParams = 2;
kern.decay = 1;
kern.variance = 1;

% Constrains parameters positive for optimisation.
kern.transforms.index = [1 2];
kern.transforms.type = optimiDefaultConstraint('positive');

% Stationarity and range of inputs
kern.isStationary = true;
kern.positiveTime = false;
