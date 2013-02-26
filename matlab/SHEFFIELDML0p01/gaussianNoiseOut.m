function y = gaussianNoiseOut(noise, mu, varsigma)

% GAUSSIANNOISEOUT Compute the output of the GAUSSIAN noise given the input mean and variance.
%
%	Description:
%
%	Y = GAUSSIANNOISEOUT(NOISE, MU, VARSIGMA) computes the ouptut for
%	the Gaussian noise given input mean and variances.
%	 Returns:
%	  Y - the output from the noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the output is computed.
%	  MU - the input mean values.
%	  VARSIGMA - the input variance values.
%	
%
%	See also
%	GAUSSIANNOISEPARAMINIT, NOISEOUT, NOISECREATE


%	Copyright (c) 2004, 2005 Neil D. Lawrence



D = size(mu, 2);
y = zeros(size(mu));
for i = 1:D
  y(:, i) = mu(:, i) + noise.bias(i);
end
