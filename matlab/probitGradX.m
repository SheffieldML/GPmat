function g = probitGradX(X, Y, model, prior)

% PROBITGRADX Gradient wrt x of log-likelihood for probit noise model.

% IVM

if size(X, 1) > 1
  error('This function only takes one data-point');
end

[mu, varsigma] = ivmPosteriorMeanVar(X, model);
[df_dx, dvarsigma_dx] = ivmPosteriorGradMeanVar(X, model);
c = Y./sqrt(1+varsigma);
u = zeros(size(c));
dlnZ_df = zeros(size(c));
u = c.*(mu+model.noise.bias);
dlnZ_df = c.*ngaussian(u)./(cummGaussian(u)); 

dlnZ_dvarsigma = -.5*c.*u.*dlnZ_df;

g = dlnZ_df*df_dx' + dlnZ_dvarsigma*dvarsigma_dx';
