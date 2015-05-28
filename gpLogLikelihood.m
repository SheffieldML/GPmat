function ll = gpLogLikelihood(model)

% GPLOGLIKELIHOOD Compute the log likelihood of a GP.
% FORMAT
% DESC computes the log likelihood of a data set given a GP model.
% ARG model : the GP model for which log likelihood is to be
% computed.
% RETURN ll : the log likelihood of the data in the GP model.
%
% SEEALSO : gpCreate, gpLogLikeGradients, modelLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2009

% GP


switch model.approx
 case 'ftc'
  % No approximation, just do a full computation on K.

  % For very high D, we use the matrix S which is M*M'
  if isfield(model, 'S')
    ll = -0.5*(model.d*model.logDetK_uu + sum(sum(model.invK_uu.* ...
                                                  model.S)));
    return;
  end
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
 case {'dtc', 'dtcvar'}
  % Deterministic training conditional
  if ~isfield(model, 'isSpherical') | model.isSpherical
    E = model.K_uf*model.m;
    EET = E*E';
    if length(model.beta)==1
      ll =  -0.5*(model.d*(-(model.N-model.k)*log(model.beta) ...
                           - model.logDetK_uu +model.logdetA) ...
                  - (sum(sum(model.Ainv.*EET)) ...
                     -sum(sum(model.m.*model.m)))*model.beta);
      if strcmp(model.approx, 'dtcvar')
        ll = ll - model.d*0.5*sum(model.diagD);
      end
    else
      error('Not implemented variable length beta yet.');
    end
  else
    ll = 0;
    for i = 1:model.d
      ind = gpDataIndices(model, i);
      e = model.K_uf(:, ind)*model.m(ind, i);
      if length(model.beta)==1
        ll = ll - 0.5*((-(model.N-model.k)*log(model.beta) ...
                        - model.logDetK_uu +model.logdetA(i)) ...
                       - (e'*model.Ainv{i}*e ...
                          -model.m(ind, i)'*model.m(ind, i))* ...
                       model.beta);
        if(isnan(ll))
          error('Log likelihood is NaN')
        end
        if strcmp(model.approx, 'dtcvar')
          error('Not implemented dtcvar for non-spherical yet.');
        end
      else
        error('Not implemented variable length beta yet.');
      end
    end
  end
 case 'fitc'
  % Fully independent training conditional.
  if ~isfield(model, 'isSpherical') | model.isSpherical
    if length(model.beta)==1
      if false
        % This is the original objective
        Dinvm = model.Dinv*model.m;
        K_ufDinvm = model.K_uf*Dinvm;
        ll = -0.5*(model.d*(sum(log(model.diagD))...
                            -(model.N-model.k)*log(model.beta) ...
                            + model.detDiff)...
                   + (sum(sum(Dinvm.*model.m))...
                      - sum(sum((model.Ainv*K_ufDinvm).*K_ufDinvm)))*model.beta);
        
        ll = ll - 0.5*model.N*model.d*log(2*pi);
      else
        % This is objective to match Ed Snelson's code
        ll =  - model.d*(sum(log(diag(model.Lm))) + 0.5*(-(model.N - model.k)*log(model.beta)+(model.N*log(2*pi)+sum(log(model.diagD)))));
        for i = 1:model.d
          ll = ll - 0.5*model.beta*(model.scaledM(:, i)'*model.scaledM(:, i) ...
                                    - model.bet(:, i)'*model.bet(:, i));
        end
      end
    else
      error('Variable length Beta not implemented yet.')
    end
  else
    if length(model.beta)==1
      if false
        ll = 0;
        for i = 1:model.d
          ind = gpDataIndices(model, i);
          Dinvm = model.Dinv{i}*model.m(ind, i);
          K_ufDinvm = model.K_uf(:, ind)*Dinvm;
          ll = ll -0.5*(sum(log(model.diagD{i})) ...
                        - (length(ind) - model.k)*log(model.beta) ...
                        + model.detDiff(i) ...
                        + (sum(sum(Dinvm.*model.m(ind, i))) ...
                           - sum(sum((model.Ainv{i}*K_ufDinvm).* ...
                                     K_ufDinvm)))*model.beta ...
                        +length(ind)*log(2*pi));
        end
      else
        % This is objective to match Ed Snelson's code
        ll = 0;
        for i = 1:model.d
          ind = gpDataIndices(model, i);
          ll =  ll - (sum(log(diag(model.Lm{i}))) ...
                      + 0.5*(-(length(ind) - model.k)*log(model.beta) ...
                             +(length(ind)*log(2*pi)+sum(log(model.diagD{i})))));
          ll = ll - 0.5*model.beta*(model.scaledM{i}'*model.scaledM{i} ...
                                    - model.bet{i}'*model.bet{i});
        end
      end
    else
      error('Variable length Beta not implemented yet.')
    end
  end
 case 'pitc'
  % Partially independent training conditional.
  if ~isfield(model, 'isSpherical') | model.isSpherical
    if length(model.beta)==1
      ll = model.d*(model.logDetA-model.logDetK_uu +model.k*log(model.beta));
      % Loop through the blocks computing each part to be added.
      K_ufDinvm = zeros(model.k, model.d);
      for i = 1:length(model.blockEnd)
        ind = gpBlockIndices(model, i);
        Dinvm{i} = model.Dinv{i}*model.m(ind, :);
        K_ufDinvm = K_ufDinvm + model.K_uf(:, ind)*Dinvm{i};
      end
      ll = ll - model.beta*sum(sum((model.Ainv*K_ufDinvm).*K_ufDinvm));
      
      for i = 1:length(model.blockEnd)
        ind = gpBlockIndices(model, i); 
        ll = ll + model.d*(model.logDetD(i) ...
                           - length(ind)*log(model.beta))...
             + model.beta*sum(sum(Dinvm{i}.*model.m(ind, :)));
      end
      ll = -0.5*ll;
      ll = ll - 0.5*model.N*model.d*log(2*pi);
    else
      error('Variable Length Beta not implemented yet.')
    end
  else
    if length(model.beta)==1
      
      ll = 0;
      for j = 1:model.d
        ll = ll + model.logDetA(j)-model.logDetK_uu + model.k*log(model.beta);
        % Loop through the blocks computing each part to be added.
        K_ufDinvm = zeros(model.k, 1);
        for i = 1:length(model.blockEnd)
          ind = gpDataIndices(model, j, i);
          Dinvm{i, j} = model.Dinv{i, j}*model.m(ind, j);
          K_ufDinvm = K_ufDinvm + model.K_uf(:, ind)*Dinvm{i, j};
        end
        ll = ll - model.beta*sum(sum((model.Ainv{i}*K_ufDinvm).*K_ufDinvm));
        
        for i = 1:length(model.blockEnd)
          ind = gpDataIndices(model, j, i);
          ll = ll + model.logDetD(i, j) ...
               - length(ind)*log(model.beta) ...
               + model.beta*sum(sum(Dinvm{i, j}.*model.m(ind, j)));
          ll = ll + length(ind)*log(2*pi);
        end
      end
      ll = -0.5*ll;
    else
      error('Variable Length Beta not implemented yet.');
    end
  end
end
if model.learnScales
  ll = ll - sum(log(model.scale));
end
ll = ll - model.d*model.N/2*log(2*pi);
