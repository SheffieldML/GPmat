function L = ngaussNoiseLikelihood(noise, mu, varsigma, y)

% NGAUSSNOISELIKELIHOOD Likelihood of the data under the NGAUSS noise model.
%
%	Description:
%
%	NGAUSSNOISELIKELIHOOD(NOISE, MU, VARSIGMA, Y) returns the likelihood
%	of a data set under the  noiseless Gaussian noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the likelihood is required.
%	  MU - input mean locations for the likelihood.
%	  VARSIGMA - input variance locations for the likelihood.
%	  Y - target locations for the likelihood.
%	
%
%	See also
%	NGAUSSNOISEPARAMINIT, NGAUSSNOISELOGLIKELIHOOD, NOISELIKELIHOOD


%	Copyright (c) 2004, 2005 Neil D. Lawrence



L = gaussianNoiseLikelihood(noise, mu, varsigma, y);