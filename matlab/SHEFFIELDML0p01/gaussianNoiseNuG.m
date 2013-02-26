function [g, nu] = gaussianNoiseNuG(noise, mu, varSigma, y)

% GAUSSIANNOISENUG Compute nu and g for GAUSSIAN noise model.
%
%	Description:
%
%	[G, NU] = GAUSSIANNOISENUG(NOISE, MU, VARSIGMA, Y, Y) computes the
%	values nu and g for the Gaussian noise given the mean and variance
%	inputs as well as the output of the noise model.
%	 Returns:
%	  G - the vector g, which is the gradient of log Z with respect to
%	   the input mean.
%	  NU - the vector nu, see equation 10 of "Extensions of the
%	   Informative Vector Machine".
%	 Arguments:
%	  NOISE - the noise structure for which the nu and g are computed.
%	  MU - input mean to the noise model.
%	  VARSIGMA - input variance to the noise model.
%	  Y - target output for the noise model.
%	  Y - target output for the noise model.
%	
%
%	See also
%	GAUSSIANNOISEPARAMINIT, NOISEUPDATENUG, NOISECREATE


%	Copyright (c) 2004, 2005 Neil D. Lawrence



D = size(y, 2);
nu = 1./(noise.sigma2+varSigma);
g = zeros(size(nu));
for i = 1:D
  g(:, i) = y(:, i) - mu(:, i) - ...
      noise.bias(i);
end
g = g.*nu;
