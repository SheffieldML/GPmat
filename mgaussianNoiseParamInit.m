function noise = mgaussianNoiseParamInit(noise, y)


% MGAUSSIANNOISEPARAMINIT MGAUSSIAN noise parameter initialisation.
% This noise model is simply a Gaussian noise that allows different
% variances when multiple outputs are used.  
% 
% SEEALSO : gaussianNoiseParamInit
%
% FORMAT
% DESC initialises the multiple output Gaussian
%  noise structure with some default parameters.
% ARG noise : the noise structure which requires initialisation.
% RETURN noise : the noise structure with the default parameters placed in.
%
% SEEALSO : noiseCreate, noiseParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


if nargin > 1
  noise.numProcess = size(y, 2);
  noise.bias = zeros(1, noise.numProcess);
  for i = 1:size(y, 2)
    noise.bias(i) = mean(y(find(~isnan(y(:, i))), i));
  end
  noise.sigma2 = repmat(1e-6, 1, noise.numProcess);
else 
  noise.bias = zeros(1, noise.numProcess);
  noise.sigma2 = repmat(1e-6, 1, noise.numProcess);
end
noise.nParams = 2*noise.numProcess;

% Can handle missing values?
noise.missing = 1;

noise.transforms.index = [noise.numProcess+1:noise.nParams];
noise.transforms.type = 'exp';
