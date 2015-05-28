function kern = wienerKernParamInit(kern)

% WIENERKERNPARAMINIT WIENER kernel parameter initialisation.
%
%	Description:
%	The Wiener covariance function is the covariance function associated with
%	the Wiener process. Increments of the Wiener process between two times are
%	sampled as if from a Gaussian distribution with mean zero and standard
%	deviation equal to the time interval.
%	
%	
%
%	KERN = WIENERKERNPARAMINIT(KERN) initialises the wiener kernel
%	structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	OUKERNCREATE, KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2009 Neil D. Lawrence


if kern.inputDimension > 1
  error('WIENER kernel only valid for one-D input.')
end

kern.variance = 1;
kern.nParams = 1;

kern.transforms.index = 1;
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = false;
kern.positiveTime = true;
