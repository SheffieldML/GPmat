function model = gpUpdateKernels(model, X, X_u, rawBeta)

% GPUPDATEKERNELS Update the kernels that are needed.

% FGPLVM

if nargin < 4
  rawBeta = [];
end
if ~isempty(rawBeta)
  fhandle = str2func([model.betaTransform 'Transform']);
  model.beta = fhandle(rawBeta, 'atox');
end
switch model.approx
 case 'ftc'
  model.K_uu = kernCompute(model.kern, X);
  [model.invK_uu, U] = pdinv(model.K_uu);
  model.logDetK_uu = logdet(model.K_uu, U);
  for i = 1:model.d
    model.innerProducts(1, i) = model.m(:, i)'*model.invK_uu...
        *model.m(:, i);
  end
  
 case {'dtc', 'fitc', 'pitc'}
  model.K_uu = kernCompute(model.kern, X_u);
  
  if ~isfield(model.kern, 'whiteVariance') | model.kern.whiteVariance == 0
    % There is no white noise term so add some jitter.
    model.K_uu = model.K_uu ...
        + sparseDiag(repmat(1e-6*mean(abs(diag(model.K_uu))), size(model.K_uu, 1), 1));
  end
  model.K_uf = kernCompute(model.kern, X_u, X);
  [model.invK_uu, U] = pdinv(model.K_uu);
  model.logDetK_uu = logdet(model.K_uu, U);
end
switch model.approx
 case 'ftc'
  % do nothing.
 case 'dtc'
  % Compute A = invBetaK_uu + K_uf*K_uf'
  K_uf2 = model.K_uf*model.K_uf';
  model.A = (1/model.beta)*model.K_uu+ K_uf2;
  % This can become unstable when K_uf2 is low rank.
  [model.Ainv, U] = pdinv(model.A);
  model.logdetA = logdet(model.A, U);
 
  % compute inner products
  for i = 1:model.d
    E = model.K_uf*model.m(:, i);    
    model.innerProducts(1, i) = ...
        model.beta*(model.m(:, i)'*model.m(:, i) ...
                    - E'*model.Ainv*E);
  end
  
 case 'fitc'
  model.diagK = kernDiagCompute(model.kern, X);
  model.diagD = (1/model.beta) + model.diagK - sum(model.K_uf.*(model.invK_uu*model.K_uf), 1)';
  model.Dinv = sparseDiag(1./model.diagD);
  K_ufDinvK_uf = model.K_uf*model.Dinv*model.K_uf';
  model.A = model.K_uu + K_ufDinvK_uf;
  % This can become unstable when K_ufDinvK_uf is low rank.
  [model.Ainv, U] = pdinv(model.A);
  model.logDetA = logdet(model.A, U);
  
  % compute inner products
  for i = 1:model.d
    Dinvm = model.Dinv*model.m(:, i);
    K_ufDinvm = model.K_uf*Dinvm;
    model.innerProducts(1, i) = Dinvm'*model.m(:, i) - K_ufDinvm'*model.Ainv*K_ufDinvm;
  end

  
 case 'pitc'
  model.A = model.K_uu;
  K_ufDinvm = zeros(model.k, model.d);
  startVal = 1;
  for i = 1:length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    model.K{i} = kernCompute(model.kern, X(ind, :));
    model.D{i} = (1/model.beta)*eye(length(ind)) + model.K{i} - ...
        model.K_uf(:, ind)'*model.invK_uu*model.K_uf(:, ind);
    [model.Dinv{i}, U] = pdinv(model.D{i});
    model.logDetD(i) = logdet(model.D{i}, U);
    K_ufDinvK_uf = model.K_uf(:, ind)*model.Dinv{i}...
        *model.K_uf(:, ind)';
    model.A = model.A + K_ufDinvK_uf;
    Dinvm{i} = model.Dinv{i}*model.m(ind, :);
    K_ufDinvm = K_ufDinvm + model.K_uf(:, ind)*Dinvm{i};
    startVal = endVal + 1;
  end
  % This can become unstable when K_ufDinvK_uf is low rank.
  [model.Ainv, U] = pdinv(model.A);
  model.logDetA = logdet(model.A, U);

  % compute inner products
  for i = 1:model.d
    model.innerProducts(1, i) = - K_ufDinvm(:, i)'*model.Ainv*K_ufDinvm(:, ...
                                                      i);
  end
  startVal = 1;
  for i = 1:length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    for j = 1:model.d
      model.innerProducts(1, j) = model.innerProducts(1, j) ...
          + Dinvm{i}(:, j)'*model.m(ind, j);
    end
    startVal = endVal + 1;
  end

 otherwise
  error('Unknown approximating criterion.')
end

