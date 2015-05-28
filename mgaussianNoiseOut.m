function y = mgaussianNoiseOut(noise, mu, varsigma)


% MGAUSSIANNOISEOUT Compute the output of the MGAUSSIAN noise given the input mean and variance.
% FORMAT
% DESC computes the ouptut for the multiple output Gaussian
% noise given input mean and variances.
% ARG noise : the noise structure for which the output is computed.
% ARG mu : the input mean values.
% ARG varSigma : the input variance values.
% RETURN y : the output from the noise model.
%
% SEEALSO : mgaussianNoiseParamInit, noiseOut, noiseCreate, 
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


D = size(mu, 2);
y = zeros(size(mu));
for i = 1:D
  y(:, i) = mu(:, i) + noise.bias(i);
end

