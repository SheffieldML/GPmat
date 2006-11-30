function [model, ll] = gpsimMapUpdateF(model, options)

% GPSIMMAPUPDATEF Update posterior mean of f.
% FORMAT
% DESC updates the stored version of the posterior mean for f in
% the GPSIMMAP model.
% ARG model : the model for which f is to be updated.
% ARG options : options vector.
% RETURN model : the model with f updated.
% RETURN ll : the log likelihood after update.
%
% SEEALSO : gpsimMapFunctionalLogLikeGradient,
% gpsimMapFunctionalLogLikeHessian, gpsimMapCreate
% 
% COPYRIGHT : Neil D. Lawrence, 2006

options = [0,  1e-4, 1e-4, 1e-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-8, 0.1, 0];

display = options(1);
tolf = options(2);
tol = options(3);
iters = options(14);
if iters == 0
  iters = 100;
end
f = gpsimMapFunctionalExtractParam(model);
llold = gpsimMapFunctionalLogLikelihood(model);

for i = 1:iters
  ll = llold - 1;
  gf = gpsimMapFunctionalLogLikeGradients(model);
  new_drn = (model.covf*gf')';
  factor = 1;
  count = 0;
  fold = f;
  while ll < llold
    f = fold + factor*new_drn;
    model = gpsimMapFunctionalExpandParam(model, f);
    ll = gpsimMapFunctionalLogLikelihood(model);
    factor = factor/2;
    count = count + 1;
    if count > 1 & display
      fprintf('gpsimMapUpdateF ... pulling back, factor %2.4f\n', factor);
    end
  end
  if ll -llold < tol & max(abs(fold - f))<tolf
    break
  end
  if display
    fprintf('Iteration %d, log likelihood %2.4f\n', i, ll);
  end
  llold = ll;
end
  