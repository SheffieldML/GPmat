function ll = gpLogLikelihood(model, Y)

% GPLOGLIKELIHOOD Compute the log likelihood of a GP.

% FGPLVM

if nargin<2
  Y = model.Y;
end

switch model.approx
 case 'ftc'
  % No approximation, just do a full computation on K.
  ll = 0;
  for i = 1:size(Y, 2)
    ll = ll -.5*model.logDetK_uu- .5*Y(:, i)'*model.invK_uu*Y(:, i);
  end
 
 case 'dtc'
  % Deterministic training conditional
  E = model.K_uf*Y;
  EET = E*E';
  ll =  -0.5*(model.d*((model.N-model.k)*log(model.sigma2) ...
                       - model.logDetK_uu +model.logdetA) ...
              - (sum(sum(model.Ainv.*EET)) ...
                 -sum(sum(Y.*Y)))/model.sigma2);
 
 case 'fitc'
  % Fully independent training conditional.
  DinvY = model.Dinv*Y;
  K_ufDinvY = model.K_uf*DinvY;
  ll = -0.5*(model.d*(sum(log(model.diagD))...
                      - model.logDetK_uu + model.logDetA) ...
             + sum(sum(DinvY.*Y)) ...
             - sum(sum((model.Ainv*K_ufDinvY).*K_ufDinvY)));
 case 'pitc'
  % Partially independent training conditional.
  ll = model.d*(model.logDetA-model.logDetK_uu);
  % Loop through the blocks computing each part to be added.
  startVal = 1;
  K_ufDinvY = zeros(model.k, model.d);
  for i = 1:length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    DinvY{i} = model.Dinv{i}*Y(ind, :);
    K_ufDinvY = K_ufDinvY + model.K_uf(:, ind)*DinvY{i};
    startVal = endVal + 1;
  end
  startVal = 1;
  ll = ll - sum(sum((model.Ainv*K_ufDinvY).*K_ufDinvY));

  for i = 1:length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    ll = ll + model.d*model.logDetD(i) ...
         + sum(sum(DinvY{i}.*Y(ind, :)));
    startVal = endVal + 1;
  end
  ll = -0.5*ll;
end
