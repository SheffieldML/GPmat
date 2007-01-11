function model = ivmEpUpdateM(model, index)

% IVMEPUPDATEM Update matrix M, L, varSigma and mu for EP.
% FORMAT
% DESC performs an EP update on the IVM model's represenations for
% a given point.
% ARG model : the mode for which the update will be done.
% ARG index : the index of the point for which the update will take
% place.
%
% SEEALSO : ivmEpUpdatePoint, ivmDowndateSites
%
% COPYRIGHT : Neil D. Lawrence, 2006

% IVM

d = length(model.I);
k = find(model.I == index);
if k == d % This point has just been included so EP update is irrelevant.
  return;
end
kvector = model.kern.Kstore(:, k);
for c = 1:length(model.Sigma)
  % Compute the kth row of  Sigma
  a = model.Sigma(c).M(:, index);
  s = kvector' - a'*model.Sigma(c).M;

  v = model.Sigma(c).Linv(k:end, k);
  sLambda_k = rocholhFactorise(v);
  
  % Update M
  model.Sigma(c).M(k:end, :) = ...
      rocholForeSub(sLambda_k, ...
		    model.Sigma(c).M(k:end, :));
  model.Sigma(c).M(k:end-1) = model.Sigma(c).M(k+1:end);

  % Update L
  t2 = model.Sigma(c).Linv(k, 1:k);
  slambda_11 = sLambda_k.u(1)*sLambda_k.v(1);
  slambda_21 = sLambda_k.v*sLambda_k.u(1);
  slambda_21(1) = [];
  sLambda_22 = sLambda_k;
  sLambda_22.s(1) = [];
  sLambda_22.u(1) = [];
  sLambda_22.v(1) = [];
  sLambda_22.n = sLambda_22.n-1;
  T3prime =  rocholForeSub(sLambda_22, model.Sigma(c).Linv(k+1:end, : ...
                                                    ));
  T3prime(:, 1:k) = T3prime(:, 1:k) ...
      - rocholForeSub(sLambda_22, slambda_21*1/slambda_11)*t2;
  v = T3prime(:, k+1:end)\T3prime(:, k);
  sVtilde = rocholFactorise(v);
  V = rocholTransMultiply(sVtilde, T3prime(:, k+1:end)')';
  T33inv = rocholTransMultiply(sLambda_22, model.Sigma(c).L(k+1:end, ...
                                                    k+1:end)')';
  invV = rocholForeSub(sVtilde, T33inv);
  model.Sigma(c).Linv(k:end-1, 1:end-1) = [T3prime(:, 1:k-1) V];
  model.Sigma(c).L(k:end-1, 1:end-1) = [-invV*T3prime(:, 1:k-1)* ...
                      model.Sigma(c).L(1:k-1, 1:k-1) invV];
  %/~
  oldVarSigma = model.varSigma;
  %~/
  model.varSigma(:, c) = model.varSigma(:, c) + ((model.nu(index, c)*s).*s)';
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

%model = ivmDowndateSites(model, index);
%model.J = [model.J index];

% Ensure nu and g now accurately reflect mu and varsigma.
model = ivmUpdateNuG(model, index);
% Update the site parameter.
model = ivmUpdateSites(model, index);
d = model.d;

% Compute the values of M for the new point
for c = 1:length(model.Sigma)
  a = model.Sigma(c).M(1:end-1, index);
  s = kvector' - a'*model.Sigma(c).M(1:end-1, ...
                                                    :);
  lValInv = sqrt(model.nu(index, c));
  %/~
  % If Nu is so low then the included data-point isn't really useful.
  if lValInv < 1e-16
    warning(['Square root of nu is ' num2str(lValInv)])
  end
  %~/
  model.Sigma(c).M(end, :) = lValInv*s;
  ainv = (-a*lValInv)'/tril(model.Sigma(c).L(1:end-1, 1:end-1));
  model.Sigma(c).L(end, :) = [a' 1/lValInv];
  model.Sigma(c).Linv(end, :) = [ainv lValInv];
  % make sure Matlab knows they are lower triangular.
  model.Sigma(c).L = tril(model.Sigma(c).L);
  model.Sigma(c).Linv = tril(model.Sigma(c).Linv);
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
% Swap the columns of the kernel matrix about.
model.kern.Kstore(:, k:end-1) = model.kern.Kstore(:, k+1:end);
model.kern.Kstore(:, end) = kvector;

% Move point to end of active set.
model.I(k:end-1) = model.I(k+1:end);
model.I(end) = index;
