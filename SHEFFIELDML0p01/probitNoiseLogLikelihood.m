function L = probitNoiseLogLikelihood(noise, mu, varsigma, y)

% PROBITNOISELOGLIKELIHOOD Log likelihood of the data under the PROBIT noise model.
%
%	Description:
%
%	PROBITNOISELOGLIKELIHOOD(NOISE, MU, VARSIGMA, Y) returns the log
%	likelihood of a data set under the  probit based classification
%	noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the log likelihood is
%	   required.
%	  MU - input mean locations for the log likelihood.
%	  VARSIGMA - input variance locations for the log likelihood.
%	  Y - target locations for the log likelihood.
%	
%
%	See also
%	PROBITNOISEPARAMINIT, PROBITNOISELIKELIHOOD, NOISELOGLIKELIHOOD


%	Copyright (c) 2004, 2005 Neil D. Lawrence



D = size(y, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end

L = sum(sum(lnCumGaussian((y.*mu)./(sqrt(noise.sigma2+varsigma)))));
