function g = gaussianGradX(X, Y, model, prior)

% GAUSSIANGRADX Gradient wrt x of log-likelihood for gaussian noise model.

% IVM


if size(X, 1) > 1
  error('This function only takes one data-point');
end

[mu, varsigma] = ivmPosteriorMeanVar(X, model);
[df_dx, dvarsigma_dx] = ivmPosteriorGradMeanVar(X, model);
nu = 1./(varsigma+model.noise.sigma2);
dlnZ_df = (Y-(mu+model.noise.bias)).*nu;

dlnZ_dvarsigma = -.5*nu+.5*dlnZ_df.*dlnZ_df;

g = dlnZ_df*df_dx' + dlnZ_dvarsigma*dvarsigma_dx';
