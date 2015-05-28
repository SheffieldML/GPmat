function y = probitNoiseOut(noise, mu, varsigma)


% PROBITNOISEOUT Compute the output of the PROBIT noise given the input mean and variance.
% FORMAT
% DESC computes the ouptut for the probit based classification
% noise given input mean and variances.
% ARG noise : the noise structure for which the output is computed.
% ARG mu : the input mean values.
% ARG varSigma : the input variance values.
% RETURN y : the output from the noise model.
%
% SEEALSO : probitNoiseParamInit, noiseOut, noiseCreate, 
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


D = size(mu, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end
y = sign(mu);
