function y = cmpndNoiseOut(noise, mu, varsigma)

% CMPNDNOISEOUT Compute the output of the CMPND noise given the input mean and variance.
%
%	Description:
%
%	Y = CMPNDNOISEOUT(NOISE, MU, VARSIGMA) computes the ouptut for the
%	compound noise given input mean and variances.
%	 Returns:
%	  Y - the output from the noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the output is computed.
%	  MU - the input mean values.
%	  VARSIGMA - the input variance values.
%	
%
%	See also
%	CMPNDNOISEPARAMINIT, NOISEOUT, NOISECREATE, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



y = zeros(size(mu));
for i = 1:length(noise.comp)
  fhandle = str2func([noise.comp{i}.type 'NoiseOut']);
  y(:, i) = fhandle(noise.comp{i},...
                    mu(:, i), ...
                    varsigma(:, i));
end

