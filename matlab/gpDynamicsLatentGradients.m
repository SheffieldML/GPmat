function gX = gpDynamicsLatentGradients(model);

% GPDYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.
% FORMAT
% DESC Compute the gradients of the log likelihood of a Gaussian 
% proces dynamics model with respect to the latent points in the
% GP-LVM model.
% ARG model : the GP model for which log likelihood is to be
% computed.
% RETURN gX : the gradients with respec to latent points.
%
% SEEALSO : gpDynamicsCreate, gpDynamicsLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2009
%
% MODIFICATIONS : Carl Henrik Ek, 2006
  
% FGPLVM

  
% The first point in each sequence cannot have dynamics thus the +
% length(model.seq)
gX = zeros(model.N+length(model.seq), model.d);
  
% First compute the regular gX for these models.
[void, gX_u, gDynX] = gpLogLikeGradients(model);

% Decide where to include gX_u.
if ~strcmp(model.approx, 'ftc') & model.fixInducing
  gDynX(model.inducingIndices, :) = gDynX(model.inducingIndices, :) + gX_u;
else
  switch model.approx
   case {'dtc', 'dtcvar', 'fitc', 'pitc'}
    if isfield(model, 'inducingPrior') & ~isempty(model.inducingPrior)
      gX_u = gX_u + priorGradient(model.inducingPrior, model.X_u);
    end
   otherwise
    % do nothing
  end
end

%/~
% switch model.approx
%  case 'ftc'
%   [void, void, gDynX] = gpLogLikeGradients(model);
%  case {'dtc', 'dtcvar', 'fitc', 'pitc'}
%   % Need to pass active set too if not using 'ftc'.
%   [void, void, gDynX] = gpLogLikeGradients(model);
% end
%~/

% Compute the indices of the input and outputs to the dynamics.
ind_in = [];
ind_out = []; 
startVal=1;
for i = 1:length(model.seq)
  endVal = model.seq(i);
  ind_in = [ind_in startVal:endVal-1];
  ind_out = [ind_out startVal+1:endVal];
  startVal = endVal + 1;
end


% Deal with the fact that X appears in the *target* for the dynamics.
switch model.approx
 case 'ftc'
  for i =1:model.q
    gX(ind_out, i) = gX(ind_out, i) - 1/model.scale(i)*model.invK_uu*model.m(:, i);
    if model.diff
      gX(ind_in, i) = gX(ind_in, i) + 1/model.scale(i)*model.invK_uu*model.m(:, i);
    end
  end
 
 case {'dtc', 'dtcvar'}
  % Deterministic training conditional.
  AinvK_uf = pdinv((1/model.beta)*model.K_uu  ...
                   + model.K_uf*model.K_uf') ...
      *model.K_uf;
  for i = 1:model.q 
    gX(ind_out, i) = gX(ind_out, i) ...
        - 1/model.scale(i)*model.beta*model.m(:, i) ...
        + 1/model.scale(i)*model.beta*model.K_uf' ...
        *(AinvK_uf*model.m(:, i));
    if model.diff
      gX(ind_in, i) = gX(ind_in, i) ...
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
    gX(ind_out, i) = gX(ind_out, i) ...
        -model.beta/model.scale(i)*sparseDiag(1./model.diagD)*model.m(:, i) ...
        + model.beta/model.scale(i)*K_ufDinv'*AinvK_ufDinv*model.m(:, i);
    if model.diff
      gX(ind_in, i) = gX(ind_in, i) ...
          + model.beta/model.scale(i)*sparseDiag(1./model.diagD)*model.m(:, i) ...
          - model.beta/model.scale(i)*K_ufDinv'*AinvK_ufDinv*model.m(:, i);
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
  i = 1;
  while i <= length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    K_ufDinv = model.K_uf(:, ind)*model.Dinv{i};
    AinvK_ufDinv = model.Ainv*K_ufDinv;
    for j = 1:model.q      
      gX(ind_out(ind), j) = gX(ind_out(ind), j) - model.beta/model.scale(j)*model.Dinv{i}*model.m(ind, j);
      if model.diff
        gX(ind_in(ind), j) = gX(ind_in(ind), j) - model.beta/model.scale(j)*model.Dinv{i}*model.m(ind, j); 
      end
      gX(ind_out, j) = gX(ind_out, j) + model.beta/model.scale(j)*gXbase'*AinvK_ufDinv*model.m(ind, j);
      if model.diff
        gX(ind_in, j) = gX(ind_in, j) - model.beta/model.scale(j)*gXbase'*AinvK_ufDinv*model.m(ind, ...
                                                          j);
      end
      
    end
    startVal = endVal + 1;
    i = i + 1;
  end
end % ends approximation cases


gX(ind_in,:) = gX(ind_in,:) + gDynX;

% If there is a prior field, use it on the first point in each sequence
if isfield(model, 'prior') &  ~isempty(model.prior)
  gX(1,:) = gX(1,:) + priorGradient(model.prior, model.X(1,:));
  for i = 2:length(model.seq)
    gX(model.seq(i-1)-1,:) = gX(model.seq(i-1)-1,:) + priorGradient(model.prior, model.X(model.seq(i-1)-i,:));
  end
end


