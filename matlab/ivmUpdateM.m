function model = ivmUpdateM(model, index)

% IVMUPDATEM Update matrix M, L, v and mu.
% FORMAT
% DESC updates the stored representations in the IVM (M, L,
% v and mu) given a new data point.
% ARG model : the model for which the represenations are to be
% updated.
% ARG index : the index of the data point that is being included.
% RETURN model : the returned model with all the representations up
% to date.
%
% SEEALSO : ivmCreate, ivmAddPoint, kernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2005
% IVM

activePoint = length(model.I)+1;

% Compute the kernel(s) at the new point.
model.kern.Kstore(:, activePoint) = kernCompute(model.kern, model.X, ...
                                                    model.X(index, :));
if isfield(model.kern, 'whiteVariance')
  model.kern.Kstore(index, activePoint) = model.kern.Kstore(index, activePoint) ...
      + model.kern.whiteVariance;
end
% Compute the values of M for the new point
for c = 1:length(model.Sigma)
  if length(model.I) > 0
    a = model.Sigma(c).M(:, index);
    s = model.kern.Kstore(:, activePoint)' - a'*model.Sigma(c).M;
    lValInv = sqrt(model.nu(index, c));
    %/~
    % If Nu is so low then the included data-point isn't really useful.
    if lValInv < 1e-16
      warning(['Square root of nu is ' num2str(lValInv)])
    end
    %~/
    model.Sigma(c).M = [model.Sigma(c).M; lValInv*s];
    ainv = (-a*lValInv)'/model.Sigma(c).L;
    model.Sigma(c).L = [[model.Sigma(c).L; a'] [zeros(length(model.I),1); 1/lValInv]];
    model.Sigma(c).Linv = [[model.Sigma(c).Linv; ainv] [zeros(length(model.I),1); lValInv]];
  else
    s = model.kern.Kstore(:, 1)';
    lValInv = sqrt(model.nu(index, c));
    model.Sigma(c).M = [lValInv*s];
    
    model.Sigma(c).L = 1/lValInv;
    model.Sigma(c).Linv = lValInv;
  end
    
  %/~
  oldVarSigma = model.varSigma(:, c);
  %~/
  model.varSigma(:, c) = model.varSigma(:, c) - ((model.nu(index, c)*s).*s)';
  %/~
  if any(model.varSigma(:, c)<0)
    minVar = min(model.varSigma(:, c));
    warning(['Minimum variance ' num2str(minVar)])
    model.varSigma(find(model.varSigma(:, c) < 0), c) = 0;
  end
  % Seems like this variance can go as low as 1e-13 in, for example, demRegression1.m
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
      minVar = min(model.varSigma(:, c));
      warning(['Minimum variance ' num2str(minVar)])
    end
    %~/
  end
end
