function model = ivmUpdateM(model, index)

% IVMUPDATEM Update matrix M, L, v and mu.

% IVM
activePoint = length(model.I)+1;
for c = 1:length(model.Sigma)
  if length(model.I) > 0
    
    model.kern.Kstore(:, activePoint) = kernCompute(model.X, model.kern, ...
                                                    model.X(index, :));
    % Add diagonal term from white noise processes
    if isfield(model.kern, 'whiteVariance')
      model.kern.Kstore(index, activePoint) = model.kern.Kstore(index, activePoint) ...
          + model.kern.whiteVariance;
    end
    if ~model.Sigma(c).robust
      a = model.Sigma(c).M(:, index);
      s = model.kern.Kstore(:, activePoint)' - a'*model.Sigma(c).M;
      lValInv = sqrt(model.nu(index, c));
      %/~
      % If Nu is so low then the included data-point isn't really useful.
      if lValInv < 1e-6
        warning(['Square root of nu is ' num2str(lValInv)])
      end
      %~/
      model.Sigma(c).M = [model.Sigma(c).M; lValInv*s];
      ainv = (-a*lValInv)'/model.Sigma(c).L;
      model.Sigma(c).L = [[model.Sigma(c).L; a'] [zeros(length(model.I),1); 1/lValInv]];
      model.Sigma(c).Linv = [[model.Sigma(c).Linv; ainv] [zeros(length(model.I),1); lValInv]];
      %/~
      %model.Sigma(c).Linv = eye(size(model.Sigma(c).L))/model.Sigma(c).L;
      %~/
    else
      a = model.Sigma(c).M(:, index);
      s = model.kern.Kstore(:, activePoint)' - a'*model.Sigma(c).M;
      sqrtBeta = sqrt(model.beta(index, c));
      a = a*sqrtBeta;
      sqrtNu = sqrt(model.nu(index, c));
      lValInv = sqrtNu/sqrtBeta;
      model.Sigma(c).M = [model.Sigma(c).M; sqrtNu*s];
      ainv = (-a*lValInv)'/model.Sigma(c).L;
      model.Sigma(c).L = [[model.Sigma(c).L; a'] [zeros(length(model.I),1); 1/lValInv]];
      model.Sigma(c).Linv = [[model.Sigma(c).Linv; ainv] [zeros(length(model.I),1); lValInv]];
    end
  else
    model.kern.Kstore(:, 1) = kernCompute(model.X, model.kern, model.X(index, :));
    
    if isfield(model.kern, 'whiteVariance')
      model.kern.Kstore(index, 1) = model.kern.Kstore(index, 1) + model.kern.whiteVariance;
    end
    if ~model.Sigma(c).robust
      s = model.kern.Kstore(:, 1)';
      lValInv = sqrt(model.nu(index, c));
      model.Sigma(c).M = [lValInv*s];
      
      model.Sigma(c).L = 1/lValInv;
      model.Sigma(c).Linv = lValInv;
    else
      s = model.kern.Kstore(:, 1)';
      sqrtNu = sqrt(model.nu(index,c));
      model.Sigma(c).M = [sqrtNu*s];
      lValInv = sqrtNu/sqrt(model.beta(index,c)); 
      % = sqrt(nu)*sqrt(beta)
      
      model.Sigma(c).L = 1/lValInv;
      model.Sigma(c).Linv = lValInv;
    end
  end
    
  %/~
  oldVarSigma = model.varSigma(:, c);
  if any(model.varSigma(:, c)<0)
    warning('Variance less than zero')
  end
  %~/
  model.varSigma(:, c) = model.varSigma(:, c) - ((model.nu(index, c)*s).*s)';
  %/~
  if any(model.varSigma(:, c)<0)
    warning('Variance less than zero')
  end
  %~/
  model.mu(:, c) = model.mu(:, c) + model.g(index, c)*s';  
end
if length(model.Sigma)==1 & size(model.y, 2)>1
  for c = 2:size(model.y, 2)
    model.varSigma(:, c) = model.varSigma(:, c) ...
        - ((model.nu(index, c)*s).*s)';
    model.mu(:, c) = model.mu(:, c) + model.g(index, c)*s'; 
    %/~
    if any(model.varSigma(:, c)<0)
      warning('Variance less than zero')
    end
    %~/
  end
end
