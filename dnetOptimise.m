function model = dnetOptimise(model, display, iters)

% DNETOPTIMISE Optimise an DNET model.
% FORMAT
% DESC optimises an mixtures of Gaussians model via the 
% expectation maximisation algorithm.
% ARG model : the model to be optimised.
% RETURN model : the optimised model.
%
% SEEALSO : mmpcaCreate, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

  
  if ~model.basisStored
    if nargin < 3
      iters = 2000;
      if nargin < 2
        display = 1;
      end
    end
    

    params = dnetExtractParam(model);
    
    options = optOptions;
    if display
      options(1) = 1;
      if length(params) <= 100
        options(9) = 1;
      end
    end
    options(14) = iters;
    
    if isfield(model, 'optimiser')
      optim = str2func(model.optimiser);
    else
      optim = str2func('scg');
    end
    
    %    if strcmp(func2str(optim), 'optimiMinimize')
    %      % Carl Rasmussen's minimize function 
    %      params = optim('dnetObjectiveGradient', params, options, model);
    %    else
    % NETLAB style optimization.
    params = optim('dnetObjective', params,  options, ...
                   'dnetGradient', model);
    %    end
    
    model = dnetExpandParam(model, params);
    model = dnetEstep(model);
  else
    diffll = 1;
    iter = 0;
    ll = dnetLowerBound(model);
    while abs(diffll)>1e-6 & iter<iters
      iter = iter + 1;
      model =  dnetEstep(model);
      if display > 1
        [ll, oldll] = boundCheck(model, ll, 'E-step');
      end
      model = dnetUpdateOutputWeights(model);
      if display > 1
        [ll, oldll] = boundCheck(model, ll, 'Weights Update');
      end
      model = dnetUpdateBeta(model);
      if display > 1
        [ll, oldll] = boundCheck(model, ll, 'Beta Update');
      end
      if display > 1
      else
        oldll = ll;
        ll = dnetLowerBound(model);
      end
      
      diffll = ll -oldll;
      if display
        fprintf('Iteration %d log-likelihood: %2.6f\n', iter, ll)
        if diffll < 0
          warning(['Log likelihood went down by ' num2str(diffll)]);
        end
      end
    end
  end
end

function [ll, oldll] = boundCheck(model, oldll, step)

% BOUNDCHECK Helper function for checking bound.

  ll = dnetLowerBound(model);
  diffll = ll - oldll;
  if ll -oldll < 0
    warning(['Log likelihood went down by ' num2str(diffll) ' in ' ...
             'step: ' step ' in dnetOptimise'])
  end
end
