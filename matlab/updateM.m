function model = updateM(model, index)

% UPDATEM Update matrix M, L, v and mu.

% IVM
activePoint = length(model.I)+1;
for c = 1:length(model.Sigma)
  if length(model.I) > 0
    
    model.kern.Kstore(:, activePoint) = kernCompute(model.X, model.kern, ...
                                                      model.X(index, :));
    if isfield(model.kern, 'whiteVariance')
      model.kern.Kstore(index, activePoint) = model.kern.Kstore(index, activePoint) ...
          + model.kern.whiteVariance;
    end
    
    a = model.Sigma(c).M(:, index);
    s = model.kern.Kstore(:, activePoint)' - a'*model.Sigma(c).M;
    sqrtNu = sqrt(model.nu(index, c));
    %/~
    % If Nu is so low then the included data-point isn't really useful.
    if sqrtNu < 1e-6
      warning(['Square root of nu is ' num2str(sqrtNu)])
    end
    %~/
    model.Sigma(c).M = [model.Sigma(c).M; sqrtNu*s];
    ainv = (-a*sqrtNu)'/model.Sigma(c).L;
    model.Sigma(c).L = [[model.Sigma(c).L; a'] [zeros(length(model.I),1); 1/sqrtNu]];
    model.Sigma(c).Linv = [[model.Sigma(c).Linv; ainv] [zeros(length(model.I),1); sqrtNu]];
    %/~
    %model.Sigma(c).Linv = eye(size(model.Sigma(c).L))/model.Sigma(c).L;
    %~/
  else
    model.kern.Kstore(:, 1) = kernCompute(model.X, model.kern, model.X(index, :));
    
    if isfield(model.kern, 'whiteVariance')
      model.kern.Kstore(index, 1) = model.kern.Kstore(index, 1) + model.kern.whiteVariance;
    end
    s = model.kern.Kstore(:, 1)';
    
    model.Sigma(c).M = [sqrt(model.nu(index, c))*s];
    
    model.Sigma(c).L = sqrt(1/model.nu(index, c));
    model.Sigma(c).Linv = 1/model.Sigma(c).L;
  end
  %/~
  oldVarSigma = model.varSigma(c, :);
  if any(model.varSigma(index, c)<0)
    warning('Variance less than zero')
  end
  %~/
  model.varSigma(:, c) = model.varSigma(:, c) - model.nu(index, c)*(s.*s)';
  %/~
  if any(model.varSigma(index, c)<0)
    warning('Variance less than zero')
  end
  %~/
  model.mu(:, c) = model.mu(:, c) + model.g(index, c)*s';
  
  %/~
  %model.varSigma(inactiveSet, c) = model.varSigma(inactiveSet, c) - model.nu(index, c)*(s(inactiveSet).*s(inactiveSet))';
  %model.mu(inactiveSet, c) = model.mu(inactiveSet, c) + model.g(index, c)*s(inactiveSet)';
  %~/
end
if length(model.Sigma)==1 & size(model.y, 2)>1
  for c = 2:size(model.y, 2)
     model.varSigma(:, c) = model.varSigma(:, c) ...
         - model.nu(index, c)*(s.*s)';
     model.mu(:, c) = model.mu(:, c) + model.g(index, c)*s'; 
  end
end