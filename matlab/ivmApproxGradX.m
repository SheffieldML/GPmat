function g = ivmApproxGradX(model, x, m, beta);

% IVMAPPROXGRADX Returns the gradient of the approxmate log-likelihood wrt x.
% FORMAT
% DESC returns the gradient of the approximate log likelihood for
% the IVM with respect to a given input position.
% ARG model : the model for which the gradient is being computed.
% ARG x : the input location for which the gradient is to be
% computed.
% ARG m : the output position where the gradient is to be computed.
% ARG beta : the output noise level for which the gradient is being
% computed.
% RETURN g : the gradient of the log likelihood with respect to x.
%
% SEEALSO : ivmPosteriorMeanVar, ivmPosteriorGradMeanVar
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% IVM

[mu, varsigma] = ivmPosteriorMeanVar(model, x);
[dmu, dvs] = ivmPosteriorGradMeanVar(model, x);


D = size(m, 2);
nu = 1./(1./beta+varsigma);
dlnZ_dmu = zeros(size(nu));
for i = 1:D
  dlnZ_dmu(:, i) = m(:, i) - mu(:, i) - model.noise.bias(i);
end
dlnZ_dmu = dlnZ_dmu.*nu;
dlnZ_dvs = 0.5*(dlnZ_dmu.*dlnZ_dmu - nu);

g = dlnZ_dmu*dmu' + dlnZ_dvs*dvs';
