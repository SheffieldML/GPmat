function model = gpUpdateAD(model, X)

% GPUPDATEAD Update the representations of A and D associated with the model.
% FORMAT
% DESC updates the representations of A and D in the model when
% called by gpUpdateKernels.
% ARG model : the model for which the representations are being
% updated.
% ARG X : the X values for which the representations are being
% computed.
% RETURN model : the model with the A and D representations
% updated.
%
% SEEALSO : gpUpdateKernels, gpExpandParam
%
% COPYRIGHT : Neil D. Lawrence 2006, 2009

% GP

if nargin < 2
  X = model.X;
end

qr_decomposition = true;

switch model.approx
 case 'ftc'
  % Compute the inner product values.
  if ~isfield(model, 'S')
    if ~isfield(model, 'isSpherical') | model.isSpherical
      for i = 1:model.d
        model.innerProducts(1, i) = model.m(:, i)'*model.invK_uu...
            *model.m(:, i);
      end
    else
      for i = 1:model.d
        ind = gpDataIndices(model, i);
        model.innerProducts(1, i) = model.m(ind, i)'*model.invK_uu{i}...
            *model.m(ind, i);
      end
    end
  end
 case {'dtc', 'dtcvar'}
  if ~isfield(model, 'isSpherical') | model.isSpherical
    % Compute A = invBetaK_uu + K_uf*K_uf'

    if ~qr_decomposition 
      % This can become unstable when K_uf2 is low rank.
      K_uf2 = model.K_uf*model.K_uf';
      model.A = (1/model.beta)*model.K_uu+ K_uf2;
      [model.Ainv, U] = pdinv(model.A);
      model.logdetA = logdet(model.A, U);
    else
      Ruu = jitChol(model.K_uu);
      sqrtBeta = sqrt(model.beta);
      [Q, R] = qr([Ruu; sqrt(model.beta)*model.K_uf'],  '0');
      RA = 1/sqrtBeta*R;
      model.A = RA'*RA;
      model.Ainv = pdinv(model.A, RA);
      model.logdetA = logdet(model.A, RA);
    end
 
    % compute inner products
    for i = 1:model.d
      if ~qr_decomposition
        E = model.K_uf*model.m(:, i);    
        model.innerProducts(1, i) = ...
            model.beta*(model.m(:, i)'*model.m(:, i) ...
                        - E'*model.Ainv*E);
      else
        Z = Q(size(model.K_uu, 1)+1:end, :)'*model.m(:, i);
        model.innerProducts(1, i) = ...
            model.beta*(model.m(:, i)'*model.m(:, i) ...
                        - Z'*Z);
      end
    end
    if strcmp(model.approx, 'dtcvar')
      model.diagD = model.beta*(model.diagK ...
        - sum(model.K_uf.*(model.invK_uu*model.K_uf), 1)');
    end
  else
    if ~model.isMissingData
      K_uf2 = model.K_uf*model.K_uf';
    end
    for i = 1:model.d
      ind = gpDataIndices(model, i);
      % Compute A = invBetaK_uu + K_uf*K_uf'
      if model.isMissingData
        K_uf2 = model.K_uf(:, ind)*model.K_uf(:, ind)';
      end
      model.A{i} = (1/model.beta)*model.K_uu+ K_uf2;
      % This can become unstable when K_uf2 is low rank.
      [model.Ainv{i}, U] = pdinv(model.A{i});
      model.logdetA(i) = logdet(model.A{i}, U);
      % compute inner products
      E = model.K_uf(:, ind)*model.m(ind, i);    
      model.innerProducts(1, i) = ...
          model.beta*(model.m(ind, i)'*model.m(ind, i) ...
                      - E'*model.Ainv{i}*E);
    end
    if strcmp(model.approx, 'dtcvar')
      error('Non spherical implementation for dtcvar not yet done.')
    end
  end
  
 case 'fitc'
  model.L = jitChol(model.K_uu)';

  if ~isfield(model, 'isSpherical') | model.isSpherical
    model.diagD = 1 + model.beta*model.diagK ...
        - model.beta*sum(model.K_uf.*(model.invK_uu*model.K_uf), 1)';
    % model.diagDdivBeta = 1/model.beta + model.diagK - sum(model.K_uf.*(model.invK_uu*model.K_uf), 1)';
    model.Dinv = sparseDiag(1./model.diagD);
    % model.betaDinv = sparseDiag(1./model.diagDdivBeta);
    K_ufDinvK_uf = model.K_uf*model.Dinv*model.K_uf';
    model.A = 1/model.beta*model.K_uu + K_ufDinvK_uf;    
    % model.betaA = model.K_uu 
    % This can become unstable when K_ufDinvK_uf is low rank.
    [model.Ainv, U] = pdinv(model.A);
    model.logDetA = logdet(model.A, U);
    % model.detDiff = model.logDetA - model.logDetK_uu;
    model.detDiff = - log(model.beta)*model.k + log(det(eye(model.k) + model.beta*K_ufDinvK_uf*model.invK_uu));
    % compute inner products
    for i = 1:model.d
      Dinvm = model.Dinv*model.m(:, i);
      K_ufDinvm = model.K_uf*Dinvm;
      model.innerProducts(1, i) = model.beta*(Dinvm'*model.m(:, i) ...
          - K_ufDinvm'*model.Ainv*K_ufDinvm);
    end
  
    % Computations from Ed's implementation.
    model.V = model.L\model.K_uf;
    model.V = model.V./repmat(sqrt(model.diagD)', model.k, 1);
    model.Am = 1/model.beta*eye(model.k)+model.V*model.V';
    model.Lm = jitChol(model.Am)';
    model.invLmV = model.Lm\model.V;
    model.scaledM = model.m./repmat(sqrt(model.diagD), 1, model.d);
    model.bet = model.invLmV*model.scaledM;
  else
    for i = 1:model.d
      ind = gpDataIndices(model, i);
      model.diagD{i} = 1 + model.beta*model.diagK(ind) ...
          - model.beta*sum(model.K_uf(:, ind).*(model.invK_uu*model.K_uf(:, ind)), 1)';
      model.Dinv{i} = sparseDiag(1./model.diagD{i});
      K_ufDinvK_uf = model.K_uf(:, ind)*model.Dinv{i}...
          *model.K_uf(:, ind)';
      model.A{i} = 1/model.beta*model.K_uu + K_ufDinvK_uf;
      % This can become unstable when K_ufDinvK_uf is low rank.
      [model.Ainv{i}, U] = pdinv(model.A{i});
      model.logDetA(i) = logdet(model.A{i}, U);
      % model.detDiff = model.logDetA - model.logDetK_uu;
      model.detDiff(i) = - log(model.beta)*model.k + log(det(eye(model.k) + model.beta*K_ufDinvK_uf*model.invK_uu));
    
      % compute inner products
      Dinvm = model.Dinv{i}*model.m(ind, i);
      K_ufDinvm = model.K_uf(:, ind)*Dinvm;
      model.innerProducts(1, i) = model.beta*(Dinvm'*model.m(ind, i) - K_ufDinvm'*model.Ainv{i}*K_ufDinvm);
    
      % Computations from Ed's implementation.
      model.V{i} = model.L\model.K_uf(:, ind);
      model.V{i} = model.V{i}./repmat(sqrt(model.diagD{i})', model.k, 1);
      model.Am{i} = 1/model.beta*eye(model.k)+model.V{i}*model.V{i}';
      model.Lm{i} = jitChol(model.Am{i})';
      model.invLmV{i} = model.Lm{i}\model.V{i};
      model.scaledM{i} = model.m(ind, i)./sqrt(model.diagD{i});
      model.bet{i} = model.invLmV{i}*model.scaledM{i};

    end
      
  end
  
 case 'pitc'
  if ~isfield(model, 'isSpherical') | model.isSpherical
    model.A = 1/model.beta*model.K_uu;
    K_ufDinvm = zeros(model.k, model.d);
    for i = 1:length(model.blockEnd)
      ind = gpBlockIndices(model, i);
      model.D{i} = eye(length(ind)) + model.beta*model.K{i} - ...
          model.beta*model.K_uf(:, ind)'*model.invK_uu*model.K_uf(:, ind);
      [model.Dinv{i}, U] = pdinv(model.D{i});
      model.logDetD(i) = logdet(model.D{i}, U);
      K_ufDinvK_uf = model.K_uf(:, ind)*model.Dinv{i}...
          *model.K_uf(:, ind)';
      model.A = model.A + K_ufDinvK_uf;
      Dinvm{i} = model.Dinv{i}*model.m(ind, :);
      K_ufDinvm = K_ufDinvm + model.K_uf(:, ind)*Dinvm{i};
    end
    % This can become unstable when K_ufDinvK_uf is low rank.
    [model.Ainv, U] = pdinv(model.A);
    model.logDetA = logdet(model.A, U);
    % compute inner products
    for i = 1:model.d
      model.innerProducts(1, i) = - model.beta*K_ufDinvm(:, i)'*model.Ainv*K_ufDinvm(:, ...
                                                        i);
    end
    for i = 1:length(model.blockEnd)
      ind = gpBlockIndices(model, i);
      for j = 1:model.d
        model.innerProducts(1, j) = model.innerProducts(1, j) ...
            + model.beta*Dinvm{i}(:, j)'*model.m(ind, j);
      end
    end
  else
    for j = 1:model.d
      model.A{j} = 1/model.beta*model.K_uu;
      K_ufDinvm = zeros(model.k, model.d);
      for i = 1:length(model.blockEnd)
        ind = gpDataIndices(model, j, i);
        model.D{i, j} = eye(length(ind)) + model.beta*model.K{i, j} - ...
            model.beta*model.K_uf(:, ind)'*model.invK_uu*model.K_uf(:, ind);
        [model.Dinv{i, j}, U] = pdinv(model.D{i, j});
        model.logDetD(i, j) = logdet(model.D{i, j}, U);
        K_ufDinvK_uf = model.K_uf(:, ind)*model.Dinv{i, j}...
            *model.K_uf(:, ind)';
        model.A{j} = model.A{j} + K_ufDinvK_uf;
        Dinvm{i}(ind, j) = model.Dinv{i, j}*model.m(ind, j);
        K_ufDinvm(:, j) = K_ufDinvm(:, j) + model.K_uf(:, ind)*Dinvm{i}(ind, ...
                                                          j);
      end
      % This can become unstable when K_ufDinvK_uf is low rank.
      [model.Ainv{j}, U] = pdinv(model.A{j});
      model.logDetA(j) = logdet(model.A{j}, U);
    end
    
    % compute inner products
    for j = 1:model.d
      model.innerProducts(1, j) = - model.beta*K_ufDinvm(:, j)'*model.Ainv{j}*K_ufDinvm(:, ...
                                                        j);
    end
    for i = 1:length(model.blockEnd)
      for j = 1:model.d
        ind = gpDataIndices(model, j, i);
        model.innerProducts(1, j) = model.innerProducts(1, j) ...
            + model.beta*Dinvm{i}(ind, j)'*model.m(ind, j);
      end
    end
    
  end
 otherwise
  error('Unknown approximating criterion.')
end

