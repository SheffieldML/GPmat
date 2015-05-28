function noise = probitNoiseParamInit(noise, y)


% PROBITNOISEPARAMINIT PROBIT noise parameter initialisation.
% The probit noise model is a classification noise model. It is based on the cumulative Gaussian. If the cumulative Gaussian is defined as
%
% \phi(x) = \int_{-\infty}^x N(z|0, 1) dz
%
% Then the probit noise model is \phi(y.*(mu+bias)/(\sqrt(sigma2 + varSigma)))
% where bias and sigma2 are parameters of the noise model, y is the true class and mu and varSigma are the input mean and variance.
%
% SEEALSO : ncnmParamInit
%
% FORMAT
% DESC initialises the probit based classification
%  noise structure with some default parameters.
% ARG noise : the noise structure which requires initialisation.
% RETURN noise : the noise structure with the default parameters placed in.
%
% SEEALSO : noiseCreate, noiseParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


if nargin > 1
  nClass1 = sum(y==1, 1);
  nClass2 = sum(y==-1, 1);
  noise.bias = invCumGaussian(nClass1./(nClass2+nClass1));
  noise.numProcess = size(y, 2);
else
  noise.bias = zeros(1, noise.numProcess);
end
noise.nParams = noise.numProcess;

% This isn't optimised, it sets the gradient of the erf.
noise.sigma2 = 1e-6;


% Can handle missing values?
noise.missing = 0;