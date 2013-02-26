function noise = gaussianNoiseParamInit(noise, y)

% GAUSSIANNOISEPARAMINIT GAUSSIAN noise parameter initialisation.
%
%	Description:
%	The Gaussian noise model is the standard noise model used for
%	regression tasks. The input mean and variance is converted to an
%	output mean and variance by first adding a bias to the mean and then
%	adding a stored variance to the input variance.
%	
%	
%
%	NOISE = GAUSSIANNOISEPARAMINIT(NOISE) initialises the Gaussian noise
%	structure with some default parameters.
%	 Returns:
%	  NOISE - the noise structure with the default parameters placed in.
%	 Arguments:
%	  NOISE - the noise structure which requires initialisation.
%	
%
%	See also
%	MGAUSSIANPARAMINIT, NOISECREATE, NOISEPARAMINIT


%	Copyright (c) 2004, 2005 Neil D. Lawrence



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