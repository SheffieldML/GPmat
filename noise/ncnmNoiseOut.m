function y = ncnmNoiseOut(noise, mu, varsigma)

% NCNMNOISEOUT Ouput from null category noise model.
% FORMAT
% DESC Gives the most likely output for the null category noise
% model for a given set of input means and variances.
% ARG noise : the noise model structure for which the output is
% calculated.
% ARG mu : the set of input means.
% ARG sigma : the set of input variances.
% RETURN y : a set of labels given the input mean and variance.
%
% SEEALSO : noiseOut, ncnmNoiseLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% NOISE

D = size(mu, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end
y = sign(mu);
