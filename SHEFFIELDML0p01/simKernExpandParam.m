function kern = simKernExpandParam(kern, params)

% SIMKERNEXPANDPARAM Create kernel structure from SIM kernel's parameters.
%
%	Description:
%
%	KERN = SIMKERNEXPANDPARAM(KERN, PARAM) returns a single input motif
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
%	See also
%	SIMKERNPARAMINIT, SIMKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2006 Neil D. Lawrence


kern.decay = params(1);
kern.inverseWidth = params(2);
if isfield(kern, 'isNegativeS') && kern.isNegativeS
  kern.sensitivity = params(3);
  kern.variance = kern.sensitivity*kern.sensitivity;
else
  kern.variance = params(3);
end
if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  kern.initialVariance = params(4);
end
