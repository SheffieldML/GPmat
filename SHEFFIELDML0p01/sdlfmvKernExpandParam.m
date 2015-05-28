function kern = sdlfmvKernExpandParam(kern, params)

% SDLFMVKERNEXPANDPARAM Pass parameters from params to SDLFMV kernel
%
%	Description:
%
%	KERN = SDLFMVKERNEXPANDPARAM(KERN, PARAM) returns a switching
%	dynamical kernel structure for velocities filled with the parameters
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
%	See also
%	SDLFMVKERNPARAMINIT, SDLFMVKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2010 Mauricio A. Alvarez


% This is just a wrapper function to call the expandParam function of the
% SDLFM structure.

kern = sdlfmKernExpandParam(kern, params);

