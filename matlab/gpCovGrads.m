function [gK_uu, gK_uf, g_Lambda, g_Lambda2] = gpCovGrads(model, Y)

% GPCOVGRADS Sparse objective function gradients wrt Covariance functions for inducing variables.

% FGPLVM

switch model.approx
 case 'dtc'
  % Deterministic training conditional.
  E = model.K_uf*Y;
  EET = E*E';
  AinvEET = model.Ainv*EET;
  AinvEETAinv = AinvEET*model.Ainv;
  gK_uu = 0.5*(model.d*(model.invK_uu-(1/model.beta)*model.Ainv) ...
          - AinvEETAinv);
  
  AinvK_uf = model.Ainv*model.K_uf;
  gK_uf = -model.d*AinvK_uf-model.beta*(AinvEET*AinvK_uf-(model.Ainv*E*Y'));
  
  g_Lambda2 = 0.5*(model.d*((model.N-model.k)/model.beta ...
                            +sum(sum(model.Ainv.*model.K_uu))/(model.beta*model.beta))...
                   +sum(sum(AinvEETAinv.*model.K_uu))/model.beta ...
                   +(trace(AinvEET)-sum(sum(Y.*Y))));
  
  fhandle = str2func([model.betaTransform 'Transform']);
  g_Lambda2 = g_Lambda2*fhandle(model.beta, 'gradfact');
  g_Lambda = [];
  
 case 'fitc'
  % Fully independent training conditonal.
  K_ufDinvK_uf = model.K_uf*model.Dinv*model.K_uf';
  E = model.K_uf*model.Dinv*Y;
  EET = E*E';
  AinvEET = model.Ainv*EET;
  diagK_fuAinvEYT = sum(model.K_uf.*(model.Ainv*E*Y'), 1)';
  AinvEETAinv = AinvEET*model.Ainv;
  diagK_ufdAinvplusAinvEETAinvK_fu = ...
      sum(model.K_uf.*((model.d*model.Ainv+AinvEETAinv)*model.K_uf), 1)';
  invK_uuK_uf = model.invK_uu*model.K_uf;
  invK_uuK_ufDinv = invK_uuK_uf*model.Dinv;
  diagYYT = sum(Y.*Y, 2);
  diagQ = -model.d*model.diagD + diagYYT ...
          + diagK_ufdAinvplusAinvEETAinvK_fu...
          -2*diagK_fuAinvEYT;
  gK_uu = 0.5*(model.d*model.invK_uu ...
               - model.d*model.Ainv - AinvEETAinv ...
               + invK_uuK_ufDinv*sparseDiag(diagQ)*invK_uuK_ufDinv');
  gK_uf = -invK_uuK_ufDinv*sparseDiag(diagQ)*model.Dinv ...      
          -model.d*model.Ainv*model.K_uf*model.Dinv ...
          -AinvEETAinv*model.K_uf*model.Dinv ...
          +model.Ainv*E*Y'*model.Dinv;
  g_Lambda = 0.5*diagQ.*1./(model.diagD.*model.diagD);
  g_Lambda2 = -sum(g_Lambda)/(model.beta*model.beta);
  fhandle = str2func([model.betaTransform 'Transform']);
  g_Lambda2 = g_Lambda2*fhandle(model.beta, 'gradfact');
 
 case 'pitc' 
  % Partially independent training conditional.
  E = zeros(model.k, model.d);
  startVal = 1;
  for i = 1:length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    E = E + model.K_uf(:, ind)*model.Dinv{i}*Y(ind, :);
    startVal = endVal + 1;
  end
  AinvE = model.Ainv*E;
  AinvEET = AinvE*E';
  AinvEETAinv = AinvEET*model.Ainv;
  startVal = 1;
  for i = 1:length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    K_fuAinvEYT = model.K_uf(:, ind)'*AinvE*Y(ind, :)';
    blockQ{i} = -model.d*model.D{i} + Y(ind, :)*Y(ind, :)' ...
        + model.K_uf(:, ind)'*(model.d*model.Ainv + AinvEETAinv)*model.K_uf(:, ind)...
        -K_fuAinvEYT - K_fuAinvEYT';
    startVal = endVal + 1;
  end
  gK_uu = model.d*model.invK_uu ...
          - model.d*model.Ainv - AinvEETAinv;
  g_Lambda2 = 0;
  gK_ufBase = -(model.d*model.Ainv + AinvEETAinv)*model.K_uf ...
      + AinvE*Y';
  
  startVal = 1;
  for i = 1:length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    invK_uuK_ufDinv = model.invK_uu*model.K_uf(:, ind)*model.Dinv{i};
    gK_uu = gK_uu + invK_uuK_ufDinv*blockQ{i}*invK_uuK_ufDinv';
    
    gK_uf(:, ind) = (gK_ufBase(:, ind) ...
                     -invK_uuK_ufDinv*blockQ{i})*model.Dinv{i};
    

    g_Lambda{i} = 0.5*model.Dinv{i}*blockQ{i}*model.Dinv{i};
    g_Lambda2 = g_Lambda2 - sum(diag((g_Lambda{i})))/(model.beta*model.beta);
    startVal = endVal + 1;
  end
  gK_uu = gK_uu*0.5;
  fhandle = str2func([model.betaTransform 'Transform']);
  g_Lambda2 = g_Lambda2*fhandle(model.beta, 'gradfact');

 otherwise
  error('Unknown approximation type');
end