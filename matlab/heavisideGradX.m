function g = heavisideGradX(X, Y, model, prior)

% HEAVISIDEGRADX Gradient wrt x of log-likelihood for heaviside noise model.

% IVM

if size(X, 1) > 1
  error('This function only takes one data-point');
end

[mu, varsigma] = ivmPosteriorMeanVar(X, model);
[df_dx, dvarsigma_dx] = ivmPosteriorGradMeanVar(X, model);
c = Y./sqrt(varsigma);
u = zeros(size(c));
dlnZ_df = zeros(size(c));
u = c.*(mu+model.noise.bias);
dlnZ_df = c.*ngaussian(u)./(cummGaussian(u) + model.noise.eta./(1-2*model.noise.eta)); 

dlnZ_dvarsigma = -.5*c.*u.*dlnZ_df;

g = dlnZ_df*df_dx' + dlnZ_dvarsigma*dvarsigma_dx';
