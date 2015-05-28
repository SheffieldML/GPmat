function model = mogOptimise(model, display, iters)

% MOGOPTIMISE Optimise an MOG model.
% FORMAT
% DESC optimises an mixtures of Gaussians model via the 
% expectation maximisation algorithm.
% ARG model : the model to be optimised.
% RETURN model : the optimised model.
%
% SEEALSO : mmpcaCreate, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

diffll = 1;
iter = 0;
ll = mogLowerBound(model);
while abs(diffll)>1e-6 & iter<iters
  iter = iter + 1;
  model =  mogEstep(model);
  if display > 1
    [ll, oldll] = boundCheck(model, ll, 'E-step');
  end
  model = mogUpdatePrior(model);
  if display > 1
    [ll, oldll] = boundCheck(model, ll, 'Prior Update');
  end
  model = mogUpdateMean(model);
  if display > 1
    [ll, oldll] = boundCheck(model, ll, 'Mean Update');
  end
  model = mogUpdateCovariance(model);
  if display > 1
    [ll, oldll] = boundCheck(model, ll, 'Covariance Update');
  end
  if display > 1
  else
    oldll = ll;
    ll = mogLowerBound(model);
  end
  diffll = ll -oldll;
  if display
    fprintf('Iteration %d log-likelihood: %2.6f\n', iter, ll)
  end
end


function [ll, oldll] = boundCheck(model, oldll, step)

% BOUNDCHECK Helper function for checking bound.

ll = mogLowerBound(model);
diffll = ll - oldll;
if ll -oldll < 0
  warning(['Log likelihood went down by ' num2str(diffll) 'in ' ...
           'step: ' step ' in mogOptimise'])
end
