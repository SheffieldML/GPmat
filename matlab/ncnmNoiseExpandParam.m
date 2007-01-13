function noise = ncnmNoiseExpandParam(noise, params)

% NCNMNOISEEXPANDPARAM Expand null category noise model's structure from param vector.
% FORMAT
% DESC returns a null category noise model structure filled with the
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
% SEEALSO : ncnmNoiseParamInit, ncnmNoiseExtractParam, noiseExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% NOISE


noise.bias = params(1:noise.numProcess);
noise.gamman = params(noise.numProcess+1);
if noise.gammaSplit
  noise.gammap = params(noise.numProcess+2);
else
  noise.gammap = noise.gamman;
end

