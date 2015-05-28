function noise = noiseExpandParam(noise, params)

% NOISEEXPANDPARAM Expand the noise model's parameters from params vector.
%
%	Description:
%
%	NOISE = NOISEEXPANDPARAM(NOISE, PARAM) returns a noise structure
%	filled with the parameters in the given vector. This is used as a
%	helper function to enable parameters to be optimised in, for
%	example, the NETLAB optimisation functions.
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
%	NOISEPARAMINIT, NOISEEXTRACTPARAM


%	Copyright (c) 2004, 2005 Neil D. Lawrence


if isfield(noise, 'transforms')
  for i = 1:length(noise.transforms)
    index = noise.transforms(i).index;
    fhandle = str2func([noise.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'atox');
  end
end
fhandle = str2func([noise.type 'NoiseExpandParam']);
noise = fhandle(noise, params);
