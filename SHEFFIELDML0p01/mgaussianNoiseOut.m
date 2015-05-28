function y = mgaussianNoiseOut(noise, mu, varsigma)

% MGAUSSIANNOISEOUT Compute the output of the MGAUSSIAN noise given the input mean and variance.
%
%	Description:
%
%	Y = MGAUSSIANNOISEOUT(NOISE, MU, VARSIGMA) computes the ouptut for
%	the multiple output Gaussian noise given input mean and variances.
%	 Returns:
%	  Y - the output from the noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the output is computed.
%	  MU - the input mean values.
%	  VARSIGMA - the input variance values.
%	
%
%	See also
%	MGAUSSIANNOISEPARAMINIT, NOISEOUT, NOISECREATE, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



D = size(mu, 2);
y = zeros(size(mu));
for i = 1:D
  y(:, i) = mu(:, i) + noise.bias(i);
end

