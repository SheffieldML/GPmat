function y = ncnmNoiseOut(noise, mu, varsigma)

% NCNMNOISEOUT Ouput from null category noise model.
%
%	Description:
%
%	Y = NCNMNOISEOUT(NOISE, MU, SIGMA) Gives the most likely output for
%	the null category noise model for a given set of input means and
%	variances.
%	 Returns:
%	  Y - a set of labels given the input mean and variance.
%	 Arguments:
%	  NOISE - the noise model structure for which the output is
%	   calculated.
%	  MU - the set of input means.
%	  SIGMA - the set of input variances.
%	
%
%	See also
%	NOISEOUT, NCNMNOISELIKELIHOOD


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


D = size(mu, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end
y = sign(mu);
