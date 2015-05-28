function L = ivmApproxLogLikelihood(model);

% IVMAPPROXLOGLIKELIHOOD Return the approximate log-likelihood for the IVM.
% FORMAT
% DESC evaluates the approximate log likelihood for the IVM. The
% approximate log likelihood involves only those points in the
% active set. It is basically the approximate likelihood for EP
% (see e.g. Kuss and Rasmussen's JMLR paper), but it is missing the
% terms that are directly dependent on the noise model.
% ARG model : the IVM model for which the approximate log
% likelihood is to be computed.
% RETURN L : the approximate log likelihood of the active set.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2004
%
% SEEALSO : ivmApproxLogLikeKernGrad, ivmApproxGradX

% IVM

x = model.X(model.I, :);
m = model.m(model.I, :);
K = kernCompute(model.kern, x);
L = 0;

if model.noise.spherical
  % there is only one value for all beta
  [invK, UC] = pdinv(K+diag(1./model.beta(model.I, 1)));
  logDetTerm = logdet(K, UC);
end
  
for i = 1:size(m, 2)
  if ~model.noise.spherical
    [invK, UC] = pdinv(K+diag(1./model.beta(model.I, i)));
    logDetTerm = logdet(K, UC);
  end
  L = L -.5*logDetTerm- .5*m(:, i)'*invK*m(:, i);
end
