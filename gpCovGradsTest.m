function model = gpCovGradsTest(model)

% GPCOVGRADSTEST Test the gradients of the likelihood wrt the covariance.
% FORMAT
% DESC tests the gradients of the covariance to ensure they are
% correct.
% ARG model : the model to be tested. 
% RETURN model : the model that was tested.
%
% SEEALSO : gpCreate, gpCovGrads
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2009

% GP

% WARNING --- this isn't testing g_Lambda in gpCovGrads
  
changeVal = 1e-6;
switch model.approx
 case 'ftc'
  
 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
  for i =1 :size(model.K_uu, 1)
    for j=1:i
      origK = model.K_uu(i, j);
      model.K_uu(i, j) = origK + changeVal;
      model.K_uu(j, i) = model.K_uu(i, j);
      [model.invK_uu, U] = pdinv(model.K_uu);
      model.logDetK_uu = logdet(model.K_uu, U);
      model = gpUpdateAD(model);
      objPlus = gpLogLikelihood(model);
      model.K_uu(i, j) = origK - changeVal;
      model.K_uu(j, i) = model.K_uu(i, j);
      [model.invK_uu, U] = pdinv(model.K_uu);
      model.logDetK_uu = logdet(model.K_uu, U);
      model = gpUpdateAD(model);
      objMinus = gpLogLikelihood(model);
      diffsK_uu(i, j) = (objPlus - objMinus)/(2*changeVal);
      diffsK_uu(j, i) = diffsK_uu(i, j);
      model.K_uu(i, j) = origK;
      model.K_uu(j, i) = origK;
      [model.invK_uu, U] = pdinv(model.K_uu);
      model.logDetK_uu = logdet(model.K_uu, U);
      model = gpUpdateAD(model);
    end
  end
  for i=1:size(model.K_uf, 1)
    for j=1:size(model.K_uf, 2)
      origK = model.K_uf(i, j);
      model.K_uf(i, j) = origK + changeVal;
      model = gpUpdateAD(model);
      objPlus = gpLogLikelihood(model);
      model.K_uf(i, j) = origK - changeVal;
      model = gpUpdateAD(model);
      objMinus = gpLogLikelihood(model);
      diffsK_uf(i, j) = (objPlus - objMinus)/(2*changeVal);
      model.K_uf(i, j) = origK;
      model = gpUpdateAD(model);
    end
  end
  
  [gK_uu, gK_uf, g_Lambda] = gpCovGrads(model, model.m);
  
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

