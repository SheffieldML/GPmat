function L = noiseExpectationLogLikelihood(noise, mu, varsigma, y);

% NOISEEXPECTATIONLOGLIKELIHOOD Return the expectation of the log likelihood.
%
%	Description:
%
%	NOISEEXPECTATIONLOGLIKELIHOOD(NOISE, MU, VARSIGMA, Y) returns the
%	expectation of the log likelihood for a gven noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the expectation of the log
%	   likelihood is required.
%	  MU - input mean locations for the likelihood.
%	  VARSIGMA - input variance locations for the likelihood.
%	  Y - target locations for the likelihood.
%	
%
%	See also
%	NOISEPARAMINIT, NOISELOGLIKELIHOOD


%	Copyright (c) 2007 Neil D. Lawrence


fhandle = str2func([noise.type 'NoiseExpectationLogLikelihood']);
L = fhandle(noise, mu, varsigma, y);


