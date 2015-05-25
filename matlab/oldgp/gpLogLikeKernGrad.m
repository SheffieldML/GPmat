function g = gpLogLikeKernGrad(model)

% GPLOGLIKEKERNGRAD Gradient of the likelihood wrt kernel parameters.

% GP

K = kernCompute(model.kern, model.X);
g = zeros(1, model.kern.nParams);

for j = 1:size(model.m, 2)
  if any(isnan(model.m(:, j)))
    tK = K;
    tX = model.X;
    tm = model.m(:, j);
    ind = find(isnan(tm));
    tX(ind, :) = [];
    tm(ind) = [];
    tK(ind, :) = [];
    tK(:, ind) = [];
    tinvK = pdinv(tK);
    covGrad = feval([model.type 'CovarianceGradient'], tinvK, tm);   
    g = g + kernGradient(model.kern, tX, covGrad);
  else
    invK = pdinv(K);
    covGrad = feval([model.type 'CovarianceGradient'], invK, model.m(:, j));
    g = g + kernGradient(model.kern, model.X, covGrad);
  end
end  
