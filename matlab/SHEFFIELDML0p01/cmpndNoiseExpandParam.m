function noise = cmpndNoiseExpandParam(noise, params)

% CMPNDNOISEEXPANDPARAM Create noise structure from CMPND noise's parameters.
%
%	Description:
%
%	NOISE = CMPNDNOISEEXPANDPARAM(NOISE, PARAM) returns a compound noise
%	structure filled with the parameters in the given vector. This is
%	used as a helper function to enable parameters to be optimised in,
%	for example, the NETLAB optimisation functions.
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
%	CMPNDNOISEPARAMINIT, CMPNDNOISEEXTRACTPARAM, NOISEEXPANDPARAM


%	Copyright (c) 2004, 2005 Neil D. Lawrence



params = params*noise.paramGroups';
startVal = 1;
endVal = 0;
for i = 1:length(noise.comp)
  endVal = endVal + noise.comp{i}.nParams;
  noise.comp{i} = noiseExpandParam(noise.comp{i}, params(1, startVal:endVal));
  startVal = endVal + 1;
end
