function L = gaussianNoiseLikelihood(noise, mu, varsigma, y)

% GAUSSIANNOISELIKELIHOOD Likelihood of the data under the GAUSSIAN noise model.
%
%	Description:
%
%	GAUSSIANNOISELIKELIHOOD(NOISE, MU, VARSIGMA, Y) returns the
%	likelihoods for data points under the  Gaussian noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the likelihood is required.
%	  MU - input mean locations for the likelihood.
%	  VARSIGMA - input variance locations for the likelihood.
%	  Y - target locations for the likelihood.
%	
%
%	See also
%	GAUSSIANNOISEPARAMINIT, GAUSSIANNOISELOGLIKELIHOOD, NOISELIKELIHOOD


%	Copyright (c) 2004, 2005 Neil D. Lawrence



N = size(y, 1);
D = size(y, 2);
varsigma = varsigma + noise.sigma2;
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end
arg = (mu - y)./sqrt(varsigma);
L = (2*pi*varsigma).^(-1/2).*exp( - .5*arg.*arg);
