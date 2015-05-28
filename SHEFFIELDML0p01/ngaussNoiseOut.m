function y = ngaussNoiseOut(noise, mu, varsigma)

% NGAUSSNOISEOUT Compute the output of the NGAUSS noise given the input mean and variance.
%
%	Description:
%
%	Y = NGAUSSNOISEOUT(NOISE, MU, VARSIGMA) computes the ouptut for the
%	noiseless Gaussian noise given input mean and variances.
%	 Returns:
%	  Y - the output from the noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the output is computed.
%	  MU - the input mean values.
%	  VARSIGMA - the input variance values.
%	
%
%	See also
%	NGAUSSNOISEPARAMINIT, NOISEOUT, NOISECREATE, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



y = gaussianNoiseOut(noise, mu, varsigma);