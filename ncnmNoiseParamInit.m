function noise = ncnmNoiseParamInit(noise, y)

% NCNMNOISEPARAMINIT null category noise model's parameter initialisation.
% The null category noise model enables semi-supervised learning
% with Gaussian processes. The approach is described in a 2004 NIPS
% paper by Lawrence and Jordan.
%
% FORMAT 
% DESC initialises the parameters of the null category noise model.
% ARG noise : the structure to initialise.
% ARG y : a set of target values.
% RETURN noise : the initialised noise structure.
%
% FORMAT 
% DESC initialises the parameters of the null category noise model.
% ARG noise : the structure to initialise.
% RETURN noise : the initialised noise structure.
%
% SEEALSO : noiseParamInit, noiseCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% NOISE

% The likelihood is not log concave.
noise.logconcave = 0;
noise.gammaSplit = 0;

if nargin > 1
  nClass1 = sum(y==1, 1);
  nClass2 = sum(y==-1, 1);
  totClass = nClass1 + nClass2;
  p1 = nClass1./totClass;
  noise.numProcess = size(y, 2);
  noise.gamman = sum(isnan(y))/length(y);
  noise.gammap = noise.gamman;
  noise.bias = invCumGaussian(p1);
else
  noise.bias = zeros(1, noise.numProcess);
  noise.gamman = 0.5;
  noise.gammap = 0.5;
end
if noise.gammaSplit
  noise.nParams = noise.numProcess+2;
else
  noise.nParams = noise.numProcess+1;
end

% Constrain noise.prior to be between 0 and 1.
if noise.gammaSplit
  noise.transforms.index = [noise.numProcess+1 noise.numProcess+2];
else
  noise.transforms.index = [noise.numProcess+1];
end
noise.transforms.type = optimiDefaultConstraint('zeroone');

% This isn't optimised, it sets the gradient of the erf.
noise.sigma2 = eps;

% Can handle missing values?
noise.missing = 1;
noise.width = 1;
