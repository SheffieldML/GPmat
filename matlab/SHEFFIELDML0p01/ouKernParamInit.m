function kern = ouKernParamInit(kern)

% OUKERNPARAMINIT Ornstein-Uhlenbeck (OU) kernel parameter initialisation.
%
%	Description:
%	The Ornstein-Uhlenbeck (OU) kernel arises driving a first-order linear
%	stochastic differential equation (SDE),
%	
%	df(t) = -theta*(f(t)-mu)dt + sigma*dw(t),
%	
%	with a standard Brownian motion or Wiener process, w(t). The covariance
%	function of the kernel is given by
%	
%	k(x_i, x_j) = 0.5*sigma2/theta * [exp(-theta*abs(x_i - x_j))
%	- exp(-theta*(x_i + x_j))]
%	
%	The parameters are sigma2, the process variance (kern.variance), and
%	theta, the inverse width of the kernel (kern.inverseWidth). The third
%	parameter of the SDE is mu, the stationary mean value of the process
%	(kern.average), which does not appear in the covariance function.
%	
%	Note that by default the covariance function obtained is non-stationary
%	due to the second term inside the brackets. However, the non-stationary
%	part can be ignored simply setting kern.isStationary = true (by default
%	kern.isStationary = false).
%	
%	There is also a multidimensional version of this kernel which includes
%	only the stationary part.
%	
%	
%
%	KERN = OUKERNPARAMINIT(KERN) initialises the Ornstein-Uhlenbeck (OU)
%	kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%	
%
%	See also
%	DEXPKERNPARAMINIT, KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2009 David Luengo
%	Copyright (c) 2009 Neil D. Lawrence

  


% Parameters of the cross-covariance function
kern.nParams = 2;
kern.decay = 1;
kern.variance = 1;

% Other parameters of the OU system
kern.average = 0;

% Constrains parameters positive for optimisation.
kern.transforms.index = [1 2];
kern.transforms.type = optimiDefaultConstraint('positive');
kern.isStationary = false;
kern.positiveTime = true;
