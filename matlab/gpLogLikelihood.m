function ll = gpLogLikelihood(model)

% GPLOGLIKELIHOOD Compute the log likelihood of a GP.
%
% ll = gpLogLikelihood(model)
%

% Copyright (c) 2006 Neil D. Lawrence
% gpLogLikelihood.m version 1.3




switch model.approx
 case 'ftc'
  % No approximation, just do a full computation on K.
  ll = 0;
  for i = 1:size(model.m, 2)
    if ~isfield(model, 'isSpherical') | model.isSpherical
      ll = ll -.5*model.logDetK_uu- .5*model.m(:, i)'*model.invK_uu*model.m(:, i);
    else
      if model.isMissingData
        m = model.m(model.indexPresent{i}, i);
      else
        m = model.m(:, i);
      end
      ll = ll - .5*model.logDetK_uu(i) - .5*m'*model.invK_uu{i}*m;
    end
  end
 case 'dtc'
  % Deterministic training conditional
  if ~isfield(model, 'isSpherical') | model.isSpherical
    E = model.K_uf*model.m;
    EET = E*E';
    ll =  -0.5*(model.d*(-(model.N-model.k)*log(model.beta) ...
                         - model.logDetK_uu +model.logdetA) ...
                - (sum(sum(model.Ainv.*EET)) ...
                   -sum(sum(model.m.*model.m)))*model.beta);
    for i = 1:size(model.m, 2)
    if model.KLCorrectionTerm
        ll = ll - 0.5.*model.beta(:,i).*sum(model.KLVariance(:,i)-model.m(:,i).*model.m(:,i));    
    end
    end
  else
    ll = 0;
    for i = 1:model.d
      ind = gpDataIndices(model, i);
      e = model.K_uf(:, ind)*model.m(ind, i);
      ll = ll - 0.5*((-(model.N-model.k)*log(model.beta) ...
                      - model.logDetK_uu +model.logdetA(i)) ...
                     - (e'*model.Ainv{i}*e ...
                   -model.m(ind, i)'*model.m(ind, i))* ...
                     model.beta);
      if model.KLCorrectionTerm
        ll = ll - 0.5.*model.beta(:,i).*sum(model.KLVariance(:,i)-model.m(:,i).*model.m(:,i));    
      end
      if(isnan(ll))
        error('Log likelihood is NaN')
      end
    end
  end
 case 'fitc'
  % Fully independent training conditional.
  if ~isfield(model, 'isSpherical') | model.isSpherical
    Dinvm = model.Dinv*model.m;
    K_ufDinvm = model.K_uf*Dinvm;
    ll = -0.5*(model.d*(sum(log(model.diagD))...
                        - model.logDetK_uu + model.logDetA) ...
               + sum(sum(Dinvm.*model.m)) ...
               - sum(sum((model.Ainv*K_ufDinvm).*K_ufDinvm)));
    for i = 1:size(model.m, 2)
    if model.KLCorrectionTerm
        ll = ll - 0.5.*model.beta.*sum(model.KLVariance(:,i)-model.m(:,i).*model.m(:,i));    
                %ll = ll - 0.5.*model.beta(:,i).*sum(model.KLVariance(:,i)-model.m(:,i).*model.m(:,i));    
    end
    end
  else
    ll = 0;
    for i = 1:model.d
      ind = gpDataIndices(model, i);
      Dinvm = model.Dinv{i}*model.m(ind, i);
      K_ufDinvm = model.K_uf(:, ind)*Dinvm;
      ll = ll -0.5*((sum(log(model.diagD{i}))...
                     - model.logDetK_uu + model.logDetA(i)) ...
                    + Dinvm'*model.m(ind, i) ...
                    - K_ufDinvm'*model.Ainv{i}*K_ufDinvm);
    end
    for i = 1:size(model.m, 2)
    if model.KLCorrectionTerm
        ll = ll - 0.5.*model.beta(:,i).*sum(model.KLVariance(:,i)-model.m(:,i).*model.m(:,i));    
    end
    end
  end
 case 'pitc'
  % Partially independent training conditional.
  if ~isfield(model, 'isSpherical') | model.isSpherical
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
    for i = 1:size(model.m, 2)
    if model.KLCorrectionTerm
        ll = ll - 0.5.*model.beta(:,i).*sum(model.KLVariance(:,i)-model.m(:,i).*model.m(:,i));    
    end
    end
  else
    ll = 0;
    for j = 1:model.d
      ll = ll + model.logDetA(j)-model.logDetK_uu;
      % Loop through the blocks computing each part to be added.
      K_ufDinvm = zeros(model.k, 1);
      for i = 1:length(model.blockEnd)
        ind = gpDataIndices(model, j, i);
        Dinvm{i, j} = model.Dinv{i, j}*model.m(ind, j);
        K_ufDinvm = K_ufDinvm + model.K_uf(:, ind)*Dinvm{i, j};
      end
      ll = ll - sum(sum((model.Ainv{i}*K_ufDinvm).*K_ufDinvm));
      
      for i = 1:length(model.blockEnd)
        ind = gpDataIndices(model, j, i);
        ll = ll + model.logDetD(i, j) ...
             + sum(sum(Dinvm{i, j}.*model.m(ind, j)));
      end
    end
    ll = -0.5*ll;
    for i = 1:size(model.m, 2)
    if model.KLCorrectionTerm
        ll = ll - 0.5.*model.beta(:,i).*sum(model.KLVariance(:,i)-model.m(:,i).*model.m(:,i));    
    end
    end
  end
   case 'nftc'
  % No approximation, just do a full computation on K.
  ll = 0;
  for i = 1:size(model.m, 2)
    if ~isfield(model, 'isSpherical') | model.isSpherical
      ll = ll -.5*model.logDetA -.5*model.m(:, i)'*model.Ainv*model.m(:, i);
      if model.KLCorrectionTerm
        ll = ll - 0.5.*sum(model.beta.*(model.KLVariance(:,i)-model.m(:,i).*model.m(:,i)));    
      end
      
    else
      if model.isMissingData
        m = model.m(model.indexPresent{i}, i);
      else
        m = model.m(:, i);
      end
      ll = ll - .5*model.logDetA(i) - .5*m'*model.Ainv{i}*m;
      if model.KLCorrectionTerm
        ll = ll - 0.5.*model.beta(:,i).*sum(model.KLVariance(:,i)-model.m(:,i).*model.m(:,i));    
      end
    end
  end
end
ll = ll - sum(log(model.scale));