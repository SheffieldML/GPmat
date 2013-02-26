function noise = ncnmNoiseExpandParam(noise, params)

% NCNMNOISEEXPANDPARAM Expand null category noise model's structure from param vector.
%
%	Description:
%
%	NOISE = NCNMNOISEEXPANDPARAM(NOISE, PARAM) returns a null category
%	noise model structure filled with the parameters in the given
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
%	NCNMNOISEPARAMINIT, NCNMNOISEEXTRACTPARAM, NOISEEXPANDPARAM


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



noise.bias = params(1:noise.numProcess);
noise.gamman = params(noise.numProcess+1);
if noise.gammaSplit
  noise.gammap = params(noise.numProcess+2);
else
  noise.gammap = noise.gamman;
end

