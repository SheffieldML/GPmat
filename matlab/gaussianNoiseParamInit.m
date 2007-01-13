function noise = gaussianNoiseParamInit(noise, y)

% GAUSSIANNOISEPARAMINIT GAUSSIAN noise parameter initialisation.
% The Gaussian noise model is the standard noise model used for
% regression tasks. The input mean and variance is converted to an
% output mean and variance by first adding a bias to the mean and then
% adding a stored variance to the input variance.
%
% SEEALSO : mgaussianParamInit
%
% FORMAT
% DESC initialises the Gaussian
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

noise.transforms.index = noise.numProcess+1;
noise.transforms.type = optimiDefaultConstraint('positive');;
noise.nParams = 1 + noise.numProcess;

% Can handle missing values?
noise.missing = 0;

% Noise model leads to constant value of beta.
noise.spherical = 1;