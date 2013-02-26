function kern = lfmKernExpandParam(kern, params)

% LFMKERNEXPANDPARAM Create kernel structure from LFM kernel's parameters.
%
%	Description:
%
%	KERN = LFMKERNEXPANDPARAM(KERN, PARAM) returns a single input motif
%	kernel structure filled with the parameters in the given vector.
%	This is used as a helper function to enable parameters to be
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
%
%	See also
%	LFMKERNPARAMINIT, LFMKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2007 Neil D. Lawrence


%	With modifications by David Luengo 2008, 2009


kern.mass = params(1);
kern.spring = params(2);
kern.damper = params(3);
kern.inverseWidth = params(4);
kern.sensitivity = params(5);

kern.alpha = kern.damper./(2*kern.mass);
kern.omega = sqrt(kern.spring./kern.mass-kern.alpha.^2);
kern.gamma = kern.alpha + j*kern.omega;

kern.zeta = kern.damper./(2*sqrt(kern.mass*kern.spring));
kern.omega_0 = sqrt(kern.spring./kern.mass);