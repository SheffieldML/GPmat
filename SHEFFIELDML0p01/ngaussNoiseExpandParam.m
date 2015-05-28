function noise = ngaussNoiseExpandParam(noise, params)

% NGAUSSNOISEEXPANDPARAM Create noise structure from NGAUSS noise's parameters.
%
%	Description:
%
%	NOISE = NGAUSSNOISEEXPANDPARAM(NOISE, PARAM) returns a noiseless
%	Gaussian noise structure filled with the parameters in the given
%	vector. This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
%	 Returns:
%	  NOISE - noise structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  NOISE - the noise structure in which the parameters are to be
%	   placed.
%	  PARAM - vector of parameters which are to be placed in the noise
%	   structure.
%	
%
%	See also
%	NGAUSSNOISEPARAMINIT, NGAUSSNOISEEXTRACTPARAM, NOISEEXPANDPARAM


%	Copyright (c) 2004, 2005 Neil D. Lawrence



noise.bias = params(1:end);