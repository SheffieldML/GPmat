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
%  
% MODIFIED : Pei Gao, 2008

% SHEFFIELDML

display = options(1);
tolf = options(2);
tol = options(3);
iters = options(14);
if iters == 0
    iters = 100;
end
f = gpsimMapFunctionalExtractParam(model);

if isfield(model, 'priorProtein') && ~isempty(model.priorProtein)
  model.consLambda = 5;
  maxLambda = 8;
  nCons = length(model.priorProtein);
else
  model.consLambda = 0;
  maxLambda = 1;
end

llold = gpsimMapFunctionalLogLikelihood(model);

for j = 1:maxLambda
  %fprintf('lambda = %d; display = %d\n', lambda, display);
  %if model.isConcave
  flag = 0;
  for i = 1:iters
    ll = llold - 1;
    gf = gpsimMapFunctionalLogLikeGradients(model);
    if model.consLambda
      consMat = zeros(size(model.invCovf));
      for k=1:nCons
        ftimeIndex = find((model.priorProteinTimes(k)-model.mapt)==0);
        consMat(ftimeIndex, ftimeIndex) = model.consLambda;
      end
      hf = inv(model.invCovf+consMat);
    else
      hf = model.covf;
    end

    new_drn = (hf*gf')';
    factor = 1;
    count = 0;
    fold = f;
    while ll < llold
        f = fold + factor*new_drn;
        model = gpsimMapFunctionalExpandParam(model, f);
        ll = gpsimMapFunctionalLogLikelihood(model);
        lldiff = ll - llold;
        count = count + 1;
        if count > 0 & display
          fprintf('gpsimMapUpdateF, lldiff: %2.4f, factor %2.4f\n', lldiff, ...
                  factor);
        end
        factor = factor/2;
        if lldiff < tol & max(abs(fold - f))<tolf
            flag = 1;
            break
        end
    end
    if display
      fprintf('Iteration %d, log likelihood %2.4f\n', i, ll);
    end
    llold = ll;
    if flag == 1
      break
    end
  end
  model.consLambda = model.consLambda*4;
end

if isempty(find(model.f))
  warning('TF is zero!');
end
%else
%  error('Not yet implemented non-concave update')
%end
