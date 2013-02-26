function L = cmpndNoiseLikelihood(noise, mu, varsigma, y)

% CMPNDNOISELIKELIHOOD Likelihood of the data under the CMPND noise model.
%
%	Description:
%
%	CMPNDNOISELIKELIHOOD(NOISE, MU, VARSIGMA, Y) returns the likelihood
%	of a data set under the  compound noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the likelihood is required.
%	  MU - input mean locations for the likelihood.
%	  VARSIGMA - input variance locations for the likelihood.
%	  Y - target locations for the likelihood.
%	
%
%	See also
%	CMPNDNOISEPARAMINIT, CMPNDNOISELOGLIKELIHOOD, NOISELIKELIHOOD


%	Copyright (c) 2004, 2005 Neil D. Lawrence



L = zeros(size(mu));
for i = 1:length(noise.comp)
  L(:, i) = noiseLikelihood(noise.comp{i},...
                mu(:, i), ...
                varsigma(:, i), ...
                y(:, i));
end