function kern = ndsimKernExpandParam(kern, params)

% NDSIMKERNEXPANDPARAM Create kernel structure from NDSIM kernel's parameters.
%
%	Description:
%
%	KERN = NDSIMKERNEXPANDPARAM(KERN, PARAM) returns a single input
%	motif kernel structure with zero decay, filled with the parameters
%	in the given vector. This is used as a helper function to enable
%	parameters to be optimised in, for example, the NETLAB optimisation
%	functions.
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
%	NDSIMKERNPARAMINIT, NDSIMKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2006 Neil D. Lawrence
%	Copyright (c) 2011 Jaakko Peltonen


kern.inverseWidth = params(1);
%fprintf(1,'NDSIM kern.inverseWidth=%f\n',kern.inverseWidth);

if isfield(kern, 'isNegativeS') && kern.isNegativeS
  kern.sensitivity = params(2);
  kern.variance = kern.sensitivity*kern.sensitivity;
else
  kern.variance = params(2);
end
if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  kern.initialVariance = params(3);
end
