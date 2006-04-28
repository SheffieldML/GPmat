function gX = gpDynamicsLatentGradients(model);

% GPDYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.
%
% gX = gpDynamicsLatentGradients(model);
%

% Copyright (c) 2006 Neil D. Lawrence
% gpDynamicsLatentGradients.m version 1.1



% gX the +1 accounts for the first point, which cannot have dynamics.
gX = zeros(model.N+1, model.d);

% First compute the regular gX for these models.
switch model.approx
 case 'ftc'
  [void, void, gDynX] = gpLogLikeGradients(...
      model);
 case {'dtc', 'fitc', 'pitc'}
  % Need to pass active set too if not using 'ftc'.
  [void, void, gDynX] = gpLogLikeGradients(...
      model);
end

% Deal with the fact that X appears in the target for the dynamics.
switch model.approx
 case 'ftc'
  for i =1:model.q
    gX(2:end, i) = gX(2:end, i) ...
        - 1/model.scale(i)*model.invK_uu*model.m(:, i);
    if model.diff
      gX(1:end-1, i) = gX(1:end-1, i) ...
          + 1/model.scale(i)*model.invK_uu*model.m(:, i);
      end
    end
   case 'dtc'
    AinvK_uf = pdinv((1/model.beta)*model.K_uu  ...
                     + model.K_uf*model.K_uf') ...
        *model.K_uf;
    for i = 1:model.q 
      gX(2:end, i) = gX(2:end, i) ...
          - 1/model.scale(i)*model.beta*model.m(:, i) ...
          + 1/model.scale(i)*model.beta*model.K_uf' ...
          *(AinvK_uf*model.m(:, i));
      if model.diff
        gX(1:end-1, i) = gX(1:end-1, i) ...
            + 1/model.scale(i)*model.beta*model.m(:, i) ...
            - 1/model.scale(i)*model.beta*model.K_uf' ...
            *(AinvK_uf*model.m(:, i));
      end
    end
    
   case 'fitc'
    % Fully independent training conditional.
    K_ufDinv = model.K_uf*model.Dinv;;
    AinvK_ufDinv = model.Ainv*K_ufDinv;
    
    for i = 1:model.q
      gX(2:end, i) = gX(2:end, i) ...
          - 1/model.scale(i)*sparseDiag(1./model.diagD)*model.m(:, i) ...
          + 1/model.scale(i)*K_ufDinv'*AinvK_ufDinv*model.m(:, i);
      if model.diff
        gX(1:end-1, i) = gX(1:end-1, i) ...
            + 1/model.scale(i)*sparseDiag(1./model.diagD)*model.m(:, i) ...
            - 1/model.scale(i)*K_ufDinv'*AinvK_ufDinv*model.m(:, i);
      end
    end
    
   case 'pitc'
    % Partially independent training conditional.
    startVal = 1;
    for i = 1:length(model.blockEnd)
      endVal = model.blockEnd(i);
      ind = startVal:endVal;
      gXbase(:, ind) =  model.K_uf(:, ind)*model.Dinv{i};
      startVal = endVal + 1;
    end
      
    startVal = 1;
    for i = 1:length(model.blockEnd)
      endVal = model.blockEnd(i);
      ind = startVal:endVal;
      K_ufDinv = model.K_uf(:, ind)*model.Dinv{i};
      AinvK_ufDinv = model.Ainv*K_ufDinv;
      for j = 1:model.q
        gX(ind+1, j) = gX(ind+1, j) ...
            - 1/model.scale(j)*model.Dinv{i}*model.m(ind, j);
        if model.diff
          gX(ind, j) = gX(ind, j) ...
              + 1/model.scale(j)*model.Dinv{i}*model.m(ind, j);
        end
        gX(2:end, j) = gX(2:end, j) + 1/model.scale(j)*gXbase'*AinvK_ufDinv*model.m(ind, j);
        if model.diff
          gX(1:end-1, j) = gX(1:end-1, j) - 1/model.scale(j)*gXbase'*AinvK_ufDinv*model.m(ind, ...
                                                            j);
        end
      end
      startVal = endVal + 1;
    end
  end
  gX(1:end-1, :) = gX(1:end-1, :) + gDynX;

  % If there is a prior field, use it on the first point.
  if isfield(model, 'prior') &  ~isempty(model.prior)
    gX(1, :) = gX(1, :) + priorGradient(model.prior, model.X(1, :)); 
  end
end