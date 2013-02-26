function y = probitNoiseOut(noise, mu, varsigma)

% PROBITNOISEOUT Compute the output of the PROBIT noise given the input mean and variance.
%
%	Description:
%
%	Y = PROBITNOISEOUT(NOISE, MU, VARSIGMA) computes the ouptut for the
%	probit based classification noise given input mean and variances.
%	 Returns:
%	  Y - the output from the noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the output is computed.
%	  MU - the input mean values.
%	  VARSIGMA - the input variance values.
%	
%
%	See also
%	PROBITNOISEPARAMINIT, NOISEOUT, NOISECREATE, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



D = size(mu, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end
y = sign(mu);
