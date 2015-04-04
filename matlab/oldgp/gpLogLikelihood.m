function L = gpLogLikelihood(model);

% GPLOGLIKELIHOOD Return the log-likelihood for a GP.

% GP

K = kernCompute(model.kern, model.X);
L = 0;

[invK, UC] = pdinv(K);

% If an entry is NaN then it is a missing value.
if ~any(isnan(model.m))  
  logDetTerm = logdet(K, UC);
  for i = 1:size(model.m, 2)
    L = L -.5*logDetTerm- .5*model.m(:, i)'*invK*model.m(:, i);
  end
else
  for i = 1:size(model.m, 2)
    
%    tInvK = invK;
%    tUC = UC;
    tK = K;
    tm = model.m(:, i);
    ind = find(isnan(tm));
%    tInvK(ind, :) = [];
%    tInvK(:, ind) = [];
%    tUC(ind, :) = [];
%    tUC(:, ind) = [];
    tK(ind, :) = [];
    tK(:, ind) = [];
    [tInvK, tUC] = pdinv(tK);
    logDetTerm = logdet(tK, tUC);
    tm(ind, :) = [];
    L = L -.5*logDetTerm- .5*tm'*tInvK*tm;
  end
end
