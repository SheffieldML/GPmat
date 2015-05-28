function g = mgaussianNoiseGradientParam(noise, mu, varsigma, y)


% MGAUSSIANNOISEGRADIENTPARAM Gradient of MGAUSSIAN noise's parameters.
% FORMAT
% DESC computes the gradient of the log Z of the multiple output Gaussian noise model with respect to the of functions with respect to the
% multiple output Gaussian
% noise's parameters. 
% ARG noise : the noise structure for which the gradients are being
% computed.
% ARG mu : the input means for which the gradients are being computed.
% ARG varSigma : the input variances for which the gradients are being computed.
% ARG y : the target values for the noise model.
% RETURN g : gradients of the log Z with respect to
% the noise parameters. The ordering of the vector should match
% that provided by the function noiseExtractParam.
%
%
% SEEALSO mgaussianNoiseParamInit, mgaussiannoiseGradVals, noiseGradientParam
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


D = size(y, 2);
u = zeros(size(y));
nu = zeros(size(y));
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
  nu(:, i) = 1./(varsigma(:, i) + noise.sigma2(i));
end

u = y - mu;

% Remove unlabelled points from the gradient.
u(find(isnan(y))) = 0;
nu(find(isnan(y)))= 0;
u = u.*nu;
gbias = sum(u, 1);
gsigma2 = -.5*sum(nu - u.*u, 1);
g = [gbias gsigma2];
