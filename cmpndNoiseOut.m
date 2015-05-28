function y = cmpndNoiseOut(noise, mu, varsigma)


% CMPNDNOISEOUT Compute the output of the CMPND noise given the input mean and variance.
% FORMAT
% DESC computes the ouptut for the compound
% noise given input mean and variances.
% ARG noise : the noise structure for which the output is computed.
% ARG mu : the input mean values.
% ARG varSigma : the input variance values.
% RETURN y : the output from the noise model.
%
% SEEALSO : cmpndNoiseParamInit, noiseOut, noiseCreate, 
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


y = zeros(size(mu));
for i = 1:length(noise.comp)
  fhandle = str2func([noise.comp{i}.type 'NoiseOut']);
  y(:, i) = fhandle(noise.comp{i},...
                    mu(:, i), ...
                    varsigma(:, i));
end

