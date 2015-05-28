function noise = ngaussNoiseParamInit(noise, y)


% NGAUSSNOISEPARAMINIT NGAUSS noise parameter initialisation.
% FORMAT
% DESC initialises the noiseless Gaussian
%  noise structure with some default parameters.
% ARG noise : the noise structure which requires initialisation.
% RETURN noise : the noise structure with the default parameters placed in.
%
% SEEALSO : noiseCreate, noiseParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


if nargin > 1
  noise.bias = mean(y);
  noise.numProcess = size(y, 2);
else 
  noise.bias = zeros(1, noise.numProcess);
end

noise.sigma2 = 1e-6;

noise.nParams = noise.numProcess;

% Can handle missing values?
noise.missing = 0;

% Noise model leads to constant value of beta.
noise.spherical = 1;