function y = orderedNoiseOut(noise, mu, varsigma)


% ORDEREDNOISEOUT Compute the output of the ORDERED noise given the input mean and variance.
% FORMAT
% DESC computes the ouptut for the ordered categorical
% noise given input mean and variances.
% ARG noise : the noise structure for which the output is computed.
% ARG mu : the input mean values.
% ARG varSigma : the input variance values.
% RETURN y : the output from the noise model.
%
% SEEALSO : orderedNoiseParamInit, noiseOut, noiseCreate, 
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


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