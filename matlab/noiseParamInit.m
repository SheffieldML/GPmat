function noise = noiseParamInit(noise, y)

% NOISEPARAMINIT Noise model's parameter initialisation.

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
