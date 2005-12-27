function model = gpUpdateKernels(model, X, X_u, rawSigma2)

% GPUPDATEKERNELS Update the kernels that are needed.

% FGPLVM

if nargin < 4
  rawSigma2 = [];
end
if ~isempty(rawSigma2)
  fhandle = str2func([model.sigma2Transform 'Transform']);
  model.sigma2 = fhandle(rawSigma2, 'atox');
end
switch model.approx
 case 'ftc'
  model.K_uu = kernCompute(model.kern, X);
  [model.invK_uu, U] = pdinv(model.K_uu);
  model.logDetK_uu = logdet(model.K_uu, U);
 
 case {'dtc', 'fitc', 'pitc'}
  model.K_uu = kernCompute(model.kern, X_u);
  
  if model.kern.whiteVariance == 0
    % Add some jitter.
    model.K_uu = model.K_uu ...
        + sparseDiag(repmat(1e-6, size(model.K_uu, 1), 1));
  end
  model.K_uf = kernCompute(model.kern, X_u, X);
  [model.invK_uu, U] = pdinv(model.K_uu);
  model.logDetK_uu = logdet(model.K_uu, U);
end
switch model.approx
 case 'ftc'
  % do nothing.
 case 'dtc'
  % Compute A = sigma2K_uu + K_uf*K_uf'
  K_uf2 = model.K_uf*model.K_uf';
  model.A = model.sigma2*model.K_uu+ K_uf2;
  [model.Ainv, U] = pdinv(model.A);
  model.logdetA = logdet(model.A, U);
 
 case 'fitc'
  model.diagK = kernDiagCompute(model.kern, X);
  model.diagD = model.sigma2 + model.diagK - sum(model.K_uf.*(model.invK_uu*model.K_uf), 1)';
  model.Dinv = sparseDiag(1./model.diagD);
  K_ufDinvK_uf = model.K_uf*model.Dinv*model.K_uf';
  model.A = model.K_uu + K_ufDinvK_uf;
  [model.Ainv, U] = pdinv(model.A);
  model.logDetA = logdet(model.A, U);
  
 case 'pitc'
  model.A = model.K_uu;
  startVal = 1;
  for i = 1:length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    model.K{i} = kernCompute(model.kern, X(ind, :));
    model.D{i} = model.sigma2*eye(length(ind)) + model.K{i} - ...
        model.K_uf(:, ind)'*model.invK_uu*model.K_uf(:, ind);
    [model.Dinv{i}, U] = pdinv(model.D{i});
    model.logDetD(i) = logdet(model.D{i}, U);
    K_ufDinvK_uf = model.K_uf(:, ind)*model.Dinv{i}...
        *model.K_uf(:, ind)';
    model.A = model.A + K_ufDinvK_uf;
    startVal = endVal + 1;
  end
  [model.Ainv, U] = pdinv(model.A);
  model.logDetA = logdet(model.A, U);
 
 otherwise
  error('Unknown approximating criterion.')
end

