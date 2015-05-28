function y = noiseOut(noise, mu, varsigma);

% NOISEOUT Give the output of the noise model given the mean and variance.
%
%	Description:
%
%	Y = NOISEOUT(NOISE, MU, VARSIGMA) computes the ouptut for the given
%	noise given input mean and variances.
%	 Returns:
%	  Y - the output from the noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the output is computed.
%	  MU - the input mean values.
%	  VARSIGMA - the input variance values.
%	
%
%	See also
%	NOISEPARAMINIT, NOISECREATE


%	Copyright (c) 2004, 2005 Neil D. Lawrence


fhandle = str2func([noise.type 'NoiseOut']);
y = fhandle(noise, mu, varsigma);
