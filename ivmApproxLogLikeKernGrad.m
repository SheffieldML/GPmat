function g = ivmApproxLogLikeKernGrad(model)

% IVMAPPROXLOGLIKEKERNGRAD Gradient of the approximate likelihood wrt kernel parameters.
% FORMAT
% DESC computes the gradient of the approximate log likelihood with
% respect to any kernel parameters.
% ARG model : the model for which the gradients are to be computed.
% RETURN g : gradient of the log likelihood with respect to the
% kernel parameters
% 
% SEEALSO : kernGradient, ivmCovarianceGradient, ivmCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% IVM

x = model.X(model.I, :);
m = model.m(model.I, :);
K = kernCompute(model.kern, x);
g = zeros(1, model.kern.nParams);

if model.noise.spherical
  % there is only one value for all beta
  invK = pdinv(K+diag(1./model.beta(model.I, 1)));
end

for j = 1:size(m, 2)
  if ~model.noise.spherical
    invK = pdinv(K+diag(1./model.beta(model.I, j)));
  end
  fhandle = str2func([model.type 'CovarianceGradient']);
  covGrad = fhandle(invK, m(:, j));
  g = g + kernGradient(model.kern, x, covGrad);
end  
