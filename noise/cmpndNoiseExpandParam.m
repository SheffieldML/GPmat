function noise = cmpndNoiseExpandParam(noise, params)


% CMPNDNOISEEXPANDPARAM Create noise structure from CMPND noise's parameters.
% FORMAT
% DESC returns a compound noise structure filled with the
% parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG noise : the noise structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% noise structure.
% RETURN noise : noise structure with the given parameters in the
% relevant locations.
%
% SEEALSO : cmpndNoiseParamInit, cmpndNoiseExtractParam, noiseExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


params = params*noise.paramGroups';
startVal = 1;
endVal = 0;
for i = 1:length(noise.comp)
  endVal = endVal + noise.comp{i}.nParams;
  noise.comp{i} = noiseExpandParam(noise.comp{i}, params(1, startVal:endVal));
  startVal = endVal + 1;
end
