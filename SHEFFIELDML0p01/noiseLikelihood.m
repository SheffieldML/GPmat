function L = noiseLikelihood(noise, mu, varsigma, y);

% NOISELIKELIHOOD Return the likelihood for each point under the noise model.
%
%	Description:
%
%	NOISELIKELIHOOD(NOISE, MU, VARSIGMA, Y) returns the likelihoods for
%	data points under the given noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the likelihood is required.
%	  MU - input mean locations for the likelihood.
%	  VARSIGMA - input variance locations for the likelihood.
%	  Y - target locations for the likelihood.
%	
%
%	See also
%	NOISEPARAMINIT, NOISELOGLIKELIHOOD


%	Copyright (c) 2004, 2005 Neil D. Lawrence


fhandle = str2func([noise.type 'NoiseLikelihood']);
L = fhandle(noise, mu, varsigma, y);
