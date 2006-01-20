function ll = gpLogLikelihood(model)

% GPLOGLIKELIHOOD Compute the log likelihood of a GP.

% FGPLVM


switch model.approx
 case 'ftc'
  % No approximation, just do a full computation on K.
  ll = 0;
  for i = 1:size(model.m, 2)
    ll = ll -.5*model.logDetK_uu- .5*model.m(:, i)'*model.invK_uu*model.m(:, i);
  end
 
 case 'dtc'
  % Deterministic training conditional
  E = model.K_uf*model.m;
  EET = E*E';
  ll =  -0.5*(model.d*(-(model.N-model.k)*log(model.beta) ...
                       - model.logDetK_uu +model.logdetA) ...
              - (sum(sum(model.Ainv.*EET)) ...
                 -sum(sum(model.m.*model.m)))*model.beta);
 
 case 'fitc'
  % Fully independent training conditional.
  Dinvm = model.Dinv*model.m;
  K_ufDinvm = model.K_uf*Dinvm;
  ll = -0.5*(model.d*(sum(log(model.diagD))...
                      - model.logDetK_uu + model.logDetA) ...
             + sum(sum(Dinvm.*model.m)) ...
             - sum(sum((model.Ainv*K_ufDinvm).*K_ufDinvm)));
 case 'pitc'
  % Partially independent training conditional.
  ll = model.d*(model.logDetA-model.logDetK_uu);
  % Loop through the blocks computing each part to be added.
  startVal = 1;
  K_ufDinvm = zeros(model.k, model.d);
  for i = 1:length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    Dinvm{i} = model.Dinv{i}*model.m(ind, :);
    K_ufDinvm = K_ufDinvm + model.K_uf(:, ind)*Dinvm{i};
    startVal = endVal + 1;
  end
  startVal = 1;
  ll = ll - sum(sum((model.Ainv*K_ufDinvm).*K_ufDinvm));

  for i = 1:length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    ll = ll + model.d*model.logDetD(i) ...
         + sum(sum(Dinvm{i}.*model.m(ind, :)));
    startVal = endVal + 1;
  end
  ll = -0.5*ll;
end
ll = ll - sum(log(model.scale));