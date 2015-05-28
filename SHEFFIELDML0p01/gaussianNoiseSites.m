function [m, beta] = gaussianNoiseSites(noise, g, nu, mu, varSigma, y)

% GAUSSIANNOISESITES Update the site parameters for the GAUSSIAN noise mode.
%
%	Description:
%
%	[M, BETA] = GAUSSIANNOISESITES(NOISE, G, NU, MU, VARSIGMA, Y)
%	updates the site parameters for the Gaussian noise model.
%	 Returns:
%	  M - the site mean parameters.
%	  BETA - the site precision parameters.
%	 Arguments:
%	  NOISE - the noise structure for which the site parameters are to
%	   be updated.
%	  G - values of g as retuned by gaussianNoiseNuG.
%	  NU - values of nu as retuned by gaussianNoiseNuG.
%	  MU - the mean value of the Gaussian input to the noise structure.
%	  VARSIGMA - the variance of the Gaussian input to the noise
%	   structure.
%	  Y - the target value.
%	
%
%	See also
%	GAUSSIANNOISEPARAMINIT, NOISEUPDATESITES


%	Copyright (c) 2004, 2005 Neil D. Lawrence


N = size(y, 1);
D = length(noise.bias);
beta = zeros(N, D);
for i = 1:size(y, 2)
  m(:, i) = y(:, i) - noise.bias(i);
end
beta = repmat(1./noise.sigma2, N, D);