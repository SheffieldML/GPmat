function g = gaussianNoiseGradientParam(noise, mu, varsigma, y)

% GAUSSIANNOISEGRADIENTPARAM Gradient of GAUSSIAN noise's parameters.
%
%	Description:
%
%	G = GAUSSIANNOISEGRADIENTPARAM(NOISE, MU, VARSIGMA, Y) computes the
%	gradient of the log Z of the Gaussian noise model with respect to
%	the of functions with respect to the Gaussian noise's parameters.
%	 Returns:
%	  G - gradients of the log probability with respect to the noise
%	   parameters. The ordering of the vector should match that provided
%	   by the function noiseExtractParam.
%	 Arguments:
%	  NOISE - the noise structure for which the gradients are being
%	   computed.
%	  MU - the input means for which the gradients are being computed.
%	  VARSIGMA - the input variances for which the gradients are being
%	   computed.
%	  Y - the target values for the noise model.
%	
%	
%
%	See also
%	% SEEALSO GAUSSIANNOISEPARAMINIT, GAUSSIANNOISEGRADVALS, NOISEGRADIENTPARAM


%	Copyright (c) 2004, 2005 Neil D. Lawrence



D = size(y, 2);
u = zeros(size(y));

for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end

u = y - mu;
nu = 1./(varsigma+noise.sigma2);
u = u.*nu;
gbias = sum(u, 1);
gsigma2 = -0.5*sum(sum(nu - u.*u, 1));
g = [gbias gsigma2];