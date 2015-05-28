function gX = gpReversibleDynamicsLatentGradients(model);

% GPREVERSIBLEDYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.

% FGPLVM

% gX the +1 accounts for the first two points, which cannot have dynamics.
gX = zeros(model.N+2, model.d);

% First compute the regular gX for these models.
switch model.approx
 case 'ftc'
  [void, void, gDynX] = gpLogLikeGradients(...
      model);
 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
  % Need to pass active set too if not using 'ftc'.
  [void, void, gDynX] = gpLogLikeGradients(...
      model);
end

% Deal with the fact that X appears in the target for the dynamics.
switch model.approx
 case 'ftc'
  for i =1:model.q/2
    gX(3:end, i) = gX(3:end, i) ...
        - 1/model.scale(i)*model.invK_uu*model.m(:, i);
    gX(2:end-1, i) = gX(2:end-1, i) ...
        + 1/model.scale(i)*model.invK_uu*model.m(:, i);
  end
 case {'dtc', 'dtcvar'}
  AinvK_uf = pdinv((1/model.beta)*model.K_uu  ...
                   + model.K_uf*model.K_uf') ...
      *model.K_uf;
  for i = 1:model.q/2 
    gX(3:end, i) = gX(3:end, i) ...
        - 1/model.scale(i)*model.beta*model.m(:, i) ...
        + 1/model.scale(i)*model.beta*model.K_uf' ...
        *(AinvK_uf*model.m(:, i));
    gX(2:end-1, i) = gX(2:end-1, i) ...
        + 1/model.scale(i)*model.beta*model.m(:, i) ...
        - 1/model.scale(i)*model.beta*model.K_uf' ...
        *(AinvK_uf*model.m(:, i));
  end
  
 case 'fitc'
  % Fully independent training conditional.
  K_ufDinv = model.K_uf*model.Dinv;;
  AinvK_ufDinv = model.Ainv*K_ufDinv;
  
  for i = 1:model.q/2
    gX(3:end, i) = gX(3:end, i) ...
        - model.beta/model.scale(i)*sparseDiag(1./model.diagD)*model.m(:, i) ...
        + model.beta/model.scale(i)*K_ufDinv'*AinvK_ufDinv*model.m(:, i);
    gX(2:end-1, i) = gX(2:end-1, i) ...
        + model.beta/model.scale(i)*sparseDiag(1./model.diagD)*model.m(:, i) ...
        - model.beta/model.scale(i)*K_ufDinv'*AinvK_ufDinv*model.m(:, i);
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
    for j = 1:model.q/2
      gX(ind+2, j) = gX(ind+2, j) ...
          - model.beta/model.scale(j)*model.Dinv{i}*model.m(ind, j);
      gX(ind+1, j) = gX(ind+1, j) ...
          + model.beta/model.scale(j)*model.Dinv{i}*model.m(ind, j);
      gX(3:end, j) = gX(3:end, j) + model.beta/model.scale(j)*gXbase'*AinvK_ufDinv*model.m(ind, j);
      gX(2:end-1, j) = gX(2:end-1, j) - model.beta/model.scale(j)*gXbase'*AinvK_ufDinv*model.m(ind, ...
                                                          j);
    end
    startVal = endVal + 1;
  end
end

gX(2:end-1, :) = gX(2:end-1, :) + gDynX(:, 1:2);
% Account for velocity terms on input.
gX(1:end-2, :) = gX(1:end-2, :) - gDynX(:, 3:4);
gX(2:end-1, :) = gX(2:end-1, :) + gDynX(:, 3:4);

% If there is a prior field, use it on the two points.
if isfield(model, 'prior') &  ~isempty(model.prior)
  gX(1:2, :) = gX(1:2, :) + priorGradient(model.prior, model.X(1:2, :)); 
end
