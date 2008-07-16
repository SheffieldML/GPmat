function y = ngaussNoiseOut(noise, mu, varsigma)


% NGAUSSNOISEOUT Compute the output of the NGAUSS noise given the input mean and variance.
% FORMAT
% DESC computes the ouptut for the noiseless Gaussian
% noise given input mean and variances.
% ARG noise : the noise structure for which the output is computed.
% ARG mu : the input mean values.
% ARG varSigma : the input variance values.
% RETURN y : the output from the noise model.
%
% SEEALSO : ngaussNoiseParamInit, noiseOut, noiseCreate, 
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


y = gaussianNoiseOut(noise, mu, varsigma);