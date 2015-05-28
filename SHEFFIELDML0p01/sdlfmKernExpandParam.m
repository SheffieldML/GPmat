function kern = sdlfmKernExpandParam(kern, params)

% SDLFMKERNEXPANDPARAM Pass parameters from params to SDLFM kernel
%
%	Description:
%
%	KERN = SDLFMKERNEXPANDPARAM(KERN, PARAM) returns a switching
%	dynamical kernel structure filled with the parameters in the given
%	vector. This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
%	 Returns:
%	  KERN - kernel structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  KERN - the kernel structure in which the parameters are to be
%	   placed.
%	  PARAM - vector of parameters which are to be placed in the kernel
%	   structure.
%	
%
%	See also
%	SDLFMKERNPARAMINIT, SDLFMKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2010 Mauricio A. Alvarez


kern.mass = params(1);
kern.spring = params(2);
kern.damper = params(3);
kern.inverseWidth = reshape(params(kern.inverseWidthIndx), kern.nlfPerInt, kern.nIntervals);
kern.switchingTimes = params(kern.switchingTimesIndx);
kern.sensitivity = reshape(params(kern.sensitivityIndx), kern.nlfPerInt, kern.nIntervals);

% Include additional quantities for information

kern.alpha = kern.damper/(2*kern.mass);
kern.omega = sqrt((4*kern.mass*kern.spring - kern.damper^2)/(4*kern.mass^2));
kern.gamma = kern.alpha + j*kern.omega;
kern.zeta = kern.damper./(2*sqrt(kern.mass*kern.spring));
kern.omega_0 = sqrt(kern.spring./kern.mass);