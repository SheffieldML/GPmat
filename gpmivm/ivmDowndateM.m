function model = ivmDowndateM(model, index)

% IVMDOWNDATEM Remove point from M, L, mu and varSigma.
% FORMAT 
% DESC removes the given data point from the IVM representations,
% in particular the matrices M and L. As well as from mu and
% varSigma. Unfortunately this operation can be numerically
% unstable: if it is done too many times these represenations must
% be recomputed (see ivmComputeLandM).
% ARG model : the model from which the point is to be removed.
% ARG index : the index of the point to remove from the model.
% 
% SEEALSO : ivmRemovePoint, ivmEpUpdatePoint, ivmComputeLandM,
% ivmUpdateM, rocholhFactorise, rocholFactorise, rocholTransMultiply,
% rocholForeSub
%
% COPYRIGHT : Neil D. Lawrence, 2005

% IVM

d = model.d;
k = find(model.I == index);
for c = 1:length(model.Sigma)
  % Compute the kth row of  Sigma
  a = model.Sigma(c).M(:, index);
  s = model.kern.Kstore(:, k)' - a'*model.Sigma(c).M;
  % Compute the three vectors that form Lambda.

  %sLambda_k = rocholFactorise(model.Sigma(c).Linv(k:end, k));
  v = model.Sigma(c).Linv(k:end, k);
%  v = v/eps;
%  sLambda_k = rocholFactorise(v);
  sLambda_k = rocholhFactorise(v);
  
  model.Sigma(c).M(k:end, :) = ...
      rocholForeSub(sLambda_k, ...
		    model.Sigma(c).M(k:end, :));

  model.Sigma(c).M(k, :) = [];
  t2 = model.Sigma(c).Linv(k, 1:k);
  if sLambda_k.n > 1;
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
    %    invV = V\eye(size(V));
    model.Sigma(c).Linv(end, :) = [];
    model.Sigma(c).Linv(:, end) = [];
    model.Sigma(c).Linv(k:end, :) = [T3prime(:, 1:k-1) V];
    model.Sigma(c).L(end, :) = [];
    model.Sigma(c).L(:, end) = [];
    model.Sigma(c).L(k:end, :) = [-invV*T3prime(:, 1:k-1)*model.Sigma(c).L(1:k-1, 1:k-1) invV];
  else
    model.Sigma(c).L(:, end) = []; 
    model.Sigma(c).L(end, :) = [];
    model.Sigma(c).Linv(:, end) = []; 
    model.Sigma(c).Linv(end, :) = [];
  end
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
model.kern.Kstore(:, k) = [];
