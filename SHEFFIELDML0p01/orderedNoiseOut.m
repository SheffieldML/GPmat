function y = orderedNoiseOut(noise, mu, varsigma)

% ORDEREDNOISEOUT Compute the output of the ORDERED noise given the input mean and variance.
%
%	Description:
%
%	Y = ORDEREDNOISEOUT(NOISE, MU, VARSIGMA) computes the ouptut for the
%	ordered categorical noise given input mean and variances.
%	 Returns:
%	  Y - the output from the noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the output is computed.
%	  MU - the input mean values.
%	  VARSIGMA - the input variance values.
%	
%
%	See also
%	ORDEREDNOISEPARAMINIT, NOISEOUT, NOISECREATE, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



D = size(mu, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end

y = zeros(size(mu));
i = noise.C-2;
index = find(mu > sum(noise.widths(1:i)));
y(index) = i+1;
for i = noise.C-3:-1:0
  index = find(mu > sum(noise.widths(1:i)) ...
               & mu <= sum(noise.widths(1:i+1)));
  y(index) = i+1;
end