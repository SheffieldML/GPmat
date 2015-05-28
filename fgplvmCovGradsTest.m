function fgplvmCovGradsTest(model)

% FGPLVMCOVGRADSTEST Test the gradients of the covariance.

% FGPLVM

switch model.approx
 case 'ftc'
  
 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
  for i =1 :size(model.K_uu, 1)
    for j=1:i
      origK = model.K_uu(i, j);
      model.K_uu(i, j) = origK + 1e-6;
      model.K_uu(j, i) = model.K_uu(i, j);
      [model.invK_uu, U] = pdinv(model.K_uu);
      model.logDetK_uu = logdet(model.K_uu, U);
      model = localUpdateAD(model);
      objPlus = fgplvmLogLikelihood(model);
      model.K_uu(i, j) = origK - 1e-6;
      model.K_uu(j, i) = model.K_uu(i, j);
      [model.invK_uu, U] = pdinv(model.K_uu);
      model.logDetK_uu = logdet(model.K_uu, U);
      model = localUpdateAD(model);
      objMinus = fgplvmLogLikelihood(model);
      diffsK_uu(i, j) = (objPlus - objMinus)/2e-6;
      diffsK_uu(j, i) = diffsK_uu(i, j);
      model.K_uu(i, j) = origK;
      model.K_uu(j, i) = origK;
      [model.invK_uu, U] = pdinv(model.K_uu);
      model.logDetK_uu = logdet(model.K_uu, U);
      model = localUpdateAD(model);
    end
  end
  for i =1 :size(model.K_uf, 1)
    for j=1:size(model.K_uf, 2)
      origK = model.K_uf(i, j);
      model.K_uf(i, j) = origK + 1e-6;
      model = localUpdateAD(model);
      objPlus = fgplvmLogLikelihood(model);
      model.K_uf(i, j) = origK - 1e-6;
      model = localUpdateAD(model);
      objMinus = fgplvmLogLikelihood(model);
      diffsK_uf(i, j) = (objPlus - objMinus)/2e-6;
      model.K_uf(i, j) = origK;
      model = localUpdateAD(model);
    end
  end
  
  [gK_uu, gK_uf, g_Lambda] = fgplvmCovGrads(model, model.m);
  
  gK_uuMaxDiff = max(max(abs(2*(gK_uu-diag(diag(gK_uu))) ...
                             + diag(diag(gK_uu)) ...
                             - diffsK_uu)));
  gK_ufMaxDiff = max(max(abs(gK_uf - diffsK_uf)));
  
  fprintf('K_uu grad max diff %2.4f\n', gK_uuMaxDiff);
  if gK_uuMaxDiff > 1e-4
    disp(2*(gK_uu-diag(diag(gK_uu))) ...
         + diag(diag(gK_uu)) ...
         - diffsK_uu);
  end
  fprintf('K_uf grad max diff %2.4f\n', gK_ufMaxDiff);
  if gK_ufMaxDiff > 1e-4
    disp(gK_uf - diffsK_uf)
  end
end


function model = localUpdateAD(model)

% LOCALUPDATED Update representation of D without kernel recomputation.

switch model.approx
 
 case {'dtc', 'dtcvar'}
  K_uf2 = model.K_uf*model.K_uf';
  model.A = (1/model.beta)*model.K_uu+ K_uf2;
  [model.Ainv, U] = pdinv(model.A);
  model.logdetA = logdet(model.A, U);
  if strcmp(model.approx, 'dtcvar')
    model.diagD = model.beta*(model.diagK ...
                              - sum(model.K_uf.*(model.invK_uu*model.K_uf), 1)');
  end
  
 case 'fitc'
  model.diagD = (1/model.beta) + model.diagK - sum(model.K_uf.*(model.invK_uu*model.K_uf), 1)';
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
    blockLength = length(ind);
    model.D{i} = (1/model.beta)*eye(blockLength) + model.K{i} - ...
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
end
