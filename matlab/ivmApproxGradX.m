function g = ivmApproxGradX(model, x, m, beta);

% IVMAPPROXGRADX Returns the gradient of the approximat log-likelihood wrt x.

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
