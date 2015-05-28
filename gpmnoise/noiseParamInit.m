function noise = noiseParamInit(noise, y)

% NOISEPARAMINIT Noise model's parameter initialisation.
% FORMAT
% DESC initialises the
%  noise structure with some default parameters.
% ARG noise : the noise structure which requires initialisation.
% RETURN noise : the noise structure with the default parameters placed in.
%
% SEEALSO : noiseCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE

% If this flag is set then the noise model leads to constant values of
% beta (e.g. a Gaussian with constant variance for each data-point &
% dimension). The default setting is set here as false.

noise.spherical = 0;
noise.logconcave = 1;

fhandle = str2func([noise.type 'NoiseParamInit']);
if nargin > 1
  noise = fhandle(noise, y);
else
  noise = fhandle(noise);
end
