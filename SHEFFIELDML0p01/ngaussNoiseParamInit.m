function noise = ngaussNoiseParamInit(noise, y)

% NGAUSSNOISEPARAMINIT NGAUSS noise parameter initialisation.
%
%	Description:
%
%	NOISE = NGAUSSNOISEPARAMINIT(NOISE) initialises the noiseless
%	Gaussian noise structure with some default parameters.
%	 Returns:
%	  NOISE - the noise structure with the default parameters placed in.
%	 Arguments:
%	  NOISE - the noise structure which requires initialisation.
%	
%
%	See also
%	NOISECREATE, NOISEPARAMINIT


%	Copyright (c) 2004, 2005 Neil D. Lawrence



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