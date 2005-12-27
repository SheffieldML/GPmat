function g = fgplvmLogLikeGradients(model)

% FGPLVMLOGLIKEGRADIENTS Compute the gradients of the EZFT sparse covariance.

% FGPLVM


[gParam, gX_u, gX] = gpLogLikeGradients(model);

% Check if Dynamics kernel is being used.
if isfield(model, 'dynamics') & ~isempty(model.dynamics)
  % Dynamics kernel is being used.
  switch model.dynamics.approx
    case 'ftc'
     [gDynParam, gDynX_u, gDynX] = gpLogLikeGradients(...
         model.dynamics, ...
         model.X(1:end-1, :),...
         model.X(2:end, :));
   case {'dtc', 'fitc', 'pitc'}
    % Need to pass active set too if not using 'ftc'.
    [gDynParam, gDynX_u, gDynX] = gpLogLikeGradients(...
        model.dynamics, ...
        model.X(1:end-1, :), ...
        model.X(2:end, :), ...
        model.X_u);
  end
  
  switch model.dynamics.approx
   case 'ftc'
    for i =1:model.q
      gX(2:end, i) = gX(2:end, i) ...
          - model.dynamics.invK_uu*model.X(2:end, i);
    end
   case 'dtc'
    AinvK_uf = pdinv(model.dynamics.sigma2*model.dynamics.K_uu  ...
                     + model.dynamics.K_uf*model.dynamics.K_uf') ...
        *model.dynamics.K_uf;
    for i = 1:model.dynamics.q 
      gX(2:end, i) = gX(2:end, i) ...
          - 1./model.dynamics.sigma2*model.X(2:end, i) ...
           + 1./model.dynamics.sigma2*model.dynamics.K_uf' ...
          *(AinvK_uf*model.X(2:end, i));
    end
    
   case 'fitc'
    % Fully independent training conditional.
    K_ufDinv = model.dynamics.K_uf*model.dynamics.Dinv;;
    AinvK_ufDinv = model.dynamics.Ainv*K_ufDinv;
    
    for i = 1:model.dynamics.q
      gX(2:end, i) = gX(2:end, i) ...
          - sparseDiag(1./model.dynamics.diagD)*model.X(2:end, i) ...
          + K_ufDinv'*AinvK_ufDinv*model.X(2:end, i);
    end
    
   case 'pitc'
    % Partially independent training conditional.
    startVal = 1;
    for i = 1:length(model.dynamics.blockEnd)
      endVal = model.dynamics.blockEnd(i);
      ind = startVal:endVal;
      gXbase(:, ind) =  model.dynamics.K_uf(:, ind)*model.dynamics.Dinv{i};
      startVal = endVal + 1;
    end

    startVal = 1;
    for i = 1:length(model.dynamics.blockEnd)
      endVal = model.dynamics.blockEnd(i);
      ind = startVal:endVal;
      K_ufDinv = model.dynamics.K_uf(:, ind)*model.dynamics.Dinv{i};
      AinvK_ufDinv = model.dynamics.Ainv*K_ufDinv;
      for j = 1:model.dynamics.q
        gX(ind+1, j) = gX(ind+1, j) ...
            - model.dynamics.Dinv{i}*model.X(ind+1, j);
        gX(2:end, j) = gX(2:end, j) + gXbase'*AinvK_ufDinv*model.X(ind+1, j);
      end
      startVal = endVal + 1;
    end
  end
  
  gX(1:end-1, :) = gX(1:end-1, :) + gDynX;
  gX_u = gX_u + gDynX_u;
  if isfield(model, 'prior') &  ~isempty(model.prior)
    gX(1, :) = gX(1, :) + priorGradient(model.prior, model.X(1, :)); 
  end
elseif isfield(model, 'prior') &  ~isempty(model.prior)
  gX = gX + priorGradient(model.prior, model.X); 
end

gParam = [gX_u(:)' gParam];

% Check for back constraints.
if isfield(model, 'back')
  g_w = modelOutputGrad(model.back, model.Y);
  g_modelParams = zeros(size(g_w, 2), 1);
  for i = 1:model.q
    g_modelParams = g_modelParams + g_w(:, :, i)'*gX(:, i);
  end
  g = [g_modelParams(:)' gParam];
else
  g = [gX(:)' gParam];
end




