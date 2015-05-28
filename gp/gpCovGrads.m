function [gK_uu, gK_uf, g_Lambda, gBeta] = gpCovGrads(model, M)

% GPCOVGRADS Sparse objective function gradients wrt Covariance functions for inducing variables.
% FORMAT
% DESC gives the gradients of the log likelihood with respect to the
% components of the sparse covariance (or the full covariance for the
% ftc case). 
% ARG model : the model for which the gradients are to be computed.
% ARG M : The training data for which the computation is to be made
% RETURN gK_uu : the gradient of the likelihood with respect to the
% elements of K_uu (or in the case of the 'ftc' criterion the
% gradients with respect to the kernel).
% RETURN gK_uf : the gradient of the likelihood with respect to the
% elements of K_uf.
% RETURN gLambda : the gradient of the likelihood with respect to
% the diagonal term in the fitc approximation and the blocks of the
% pitc approximation.
% RETURN gBeta : the gradient with respect to the beta term in the
% covariance structure.
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2009
%
% SEEALSO : gpCreate, gpLogLikeGradient

% GP
qr_decomposition = false; % not yet implemented
switch model.approx
 case {'dtc', 'dtcvar'}
  % Deterministic training conditional.
  if strcmp(model.approx, 'dtcvar')
    dtcvar = true;
  else
    dtcvar = false;
  end
  if ~isfield(model, 'isSpherical') | model.isSpherical
    if ~qr_decomposition
      E = model.K_uf*M;
      EET = E*E';
      AinvEET = model.Ainv*EET;
      AinvEETAinv = AinvEET*model.Ainv;
      gK_uu = 0.5*(model.d*(model.invK_uu-(1/model.beta)*model.Ainv) ...
                   - AinvEETAinv);
      if dtcvar
        K_uuInvK_uf = model.invK_uu*model.K_uf;
        gK_uu = gK_uu - 0.5*model.d*model.beta...
                *K_uuInvK_uf*K_uuInvK_uf';
      end
      AinvK_uf = model.Ainv*model.K_uf;
      gK_uf = -model.d*AinvK_uf-model.beta*(AinvEET*AinvK_uf-(model.Ainv*E*M'));
      if dtcvar
        gK_uf = gK_uf + model.d*model.beta*K_uuInvK_uf;
      end
      gBeta = 0.5*(model.d*((model.N-model.k)/model.beta ...
                              +sum(sum(model.Ainv.*model.K_uu))/(model.beta*model.beta))...
                   +sum(sum(AinvEETAinv.*model.K_uu))/model.beta ...
                   +(trace(AinvEET)-sum(sum(M.*M))));
      if dtcvar
        gBeta = gBeta -0.5*model.d*sum(model.diagD)/model.beta;
      end
    end
    fhandle = str2func([model.betaTransform 'Transform']);
    gBeta = gBeta*fhandle(model.beta, 'gradfact');
    if dtcvar
      g_Lambda = repmat(-0.5*model.beta*model.d, 1, model.N);
    else
      g_Lambda = [];
    end
  else
    gK_uu = zeros(model.k, model.k);
    gK_uf = zeros(model.k, model.N);
    gBeta = 0;
    for i = 1:model.d
      ind = gpDataIndices(model, i);
      e = model.K_uf(:, ind)*M(ind, i);
      Ainve = model.Ainv{i}*e;
      AinveeT = Ainve*e';      
      AinveeTAinv = Ainve*Ainve';
      gK_uu = gK_uu+0.5*((model.invK_uu-(1/model.beta)*model.Ainv{i}) ...
                         - AinveeTAinv);
      
      AinvK_uf = model.Ainv{i}*model.K_uf(:, ind);
      gK_uf(:, ind) = gK_uf(:, ind) - AinvK_uf...
          -model.beta*(AinveeT*AinvK_uf-(Ainve*M(ind, i)'));
      
      gBeta = gBeta ...
          + 0.5*(((model.N-model.k)/model.beta ...
                  +sum(sum(model.Ainv{i}.*model.K_uu))/(model.beta*model.beta))...
                 +sum(sum(AinveeTAinv.*model.K_uu))/model.beta ...
                 +(trace(AinveeT)-sum(sum(M(ind, i).*M(ind, i)))));
    end
    fhandle = str2func([model.betaTransform 'Transform']);
    gBeta = gBeta*fhandle(model.beta, 'gradfact');
    g_Lambda = [];
  end
 case 'fitc'
  % Fully independent training conditonal.
  if ~isfield(model, 'isSpherical') | model.isSpherical
    E = model.K_uf*model.Dinv*M;
    EET = E*E';
    AinvE = model.Ainv*E;
    %AinvEET = model.Ainv*EET;
    diagK_fuAinvEMT = sum(model.K_uf.*(model.Ainv*E*M'), 1)';
    AinvEETAinv = AinvE*AinvE';
    diagK_ufdAinvplusAinvEETAinvK_fu = ...
        sum(model.K_uf.*((model.d*model.Ainv+model.beta*AinvEETAinv)*model.K_uf), 1)';
    invK_uuK_uf = model.invK_uu*model.K_uf;
    if true
      invK_uuK_ufDinv = invK_uuK_uf*model.Dinv;
    else
      invK_uuK_ufDinv = model.L'\model.V;
    end
    diagMMT = sum(M.*M, 2);
    diagQ = -model.d*model.diagD + model.beta*diagMMT ...
            + diagK_ufdAinvplusAinvEETAinvK_fu...
            -2*model.beta*diagK_fuAinvEMT;
    gK_uu = 0.5*(model.d*(model.invK_uu ...
                 -model.Ainv/model.beta) - AinvEETAinv ...
                 + model.beta*invK_uuK_ufDinv*sparseDiag(diagQ)*invK_uuK_ufDinv');
    gK_uf = -model.beta*invK_uuK_ufDinv*sparseDiag(diagQ)*model.Dinv ...      
            -model.d*model.Ainv*model.K_uf*model.Dinv ...
            -model.beta*AinvEETAinv*model.K_uf*model.Dinv ...
            +model.beta*model.Ainv*E*M'*model.Dinv;
    g_Lambda = (0.5*diagQ*model.beta)./(model.diagD.*model.diagD);
    gBeta = -sum(g_Lambda)/(model.beta*model.beta);
    fhandle = str2func([model.betaTransform 'Transform']);
    gBeta = gBeta*fhandle(model.beta, 'gradfact');
  else
    gK_uu = zeros(model.k, model.k);
    gK_uf = zeros(model.k, model.N);
    g_Lambda = zeros(model.N, 1);
    gBeta = 0;
    for i = 1:model.d
      ind = gpDataIndices(model, i);
      K_ufDinvK_uf = model.K_uf(:, ind)*model.Dinv{i}*model.K_uf(:, ...
                                                        ind)';
      e = model.K_uf(:, ind)*model.Dinv{i}*M(ind, i);
      Ainve = model.Ainv{i}*e;
      AinveeTAinv = Ainve*Ainve';
      diagK_fuAinveyT = sum(model.K_uf(:, ind).*(Ainve*M(ind,i)'), 1)';
      diagK_ufdAinvplusAinveeTAinvK_fu = ...
          sum(model.K_uf(:, ind).*((model.Ainv{i}+model.beta*AinveeTAinv)*model.K_uf(:, ...
                                                        ind)), 1)';
      invK_uuK_uf = model.invK_uu*model.K_uf(:, ind);
      invK_uuK_ufDinv = invK_uuK_uf*model.Dinv{i};
      diagyyT = M(ind, i).*M(ind, i);
      diagQ = -model.diagD{i} + model.beta*diagyyT ...
              + diagK_ufdAinvplusAinveeTAinvK_fu...
              -2*model.beta*diagK_fuAinveyT;
      gK_uu = gK_uu ...
              +0.5*(model.invK_uu ...
                    - model.Ainv{i}/model.beta - AinveeTAinv ...
                    + model.beta*invK_uuK_ufDinv*sparseDiag(diagQ)*invK_uuK_ufDinv');
      gK_uf(:, ind) = gK_uf(:, ind) ...
              -model.beta*invK_uuK_ufDinv*sparseDiag(diagQ)*model.Dinv{i} ...      
              -model.Ainv{i}*model.K_uf(:, ind)*model.Dinv{i} ...
              -model.beta*AinveeTAinv*model.K_uf(:, ind)*model.Dinv{i} ...
              +model.beta*Ainve*M(ind, i)'*model.Dinv{i};
      g_Lambda(ind) = g_Lambda(ind) ...
          + 0.5*model.beta*diagQ./(model.diagD{i}.*model.diagD{i});
    end
    gBeta = gBeta - sum(g_Lambda)/(model.beta*model.beta);
    fhandle = str2func([model.betaTransform 'Transform']);
    gBeta = gBeta*fhandle(model.beta, 'gradfact');
  end
 
 case 'pitc' 
  % Partially independent training conditional.
  if ~isfield(model, 'isSpherical') | model.isSpherical
    E = zeros(model.k, model.d);
    for i = 1:length(model.blockEnd)
      ind = gpBlockIndices(model, i);
      E = E + model.K_uf(:, ind)*model.Dinv{i}*M(ind, :);
    end
    AinvE = model.Ainv*E;
    AinvEET = AinvE*E';
    AinvEETAinv = AinvEET*model.Ainv;
    for i = 1:length(model.blockEnd)
      ind = gpBlockIndices(model, i);
      K_fuAinvEMT = model.beta*model.K_uf(:, ind)'*AinvE*M(ind, :)';
      blockQ{i} = -model.d*model.D{i} + model.beta*M(ind, :)*M(ind, :)' ...
          + model.K_uf(:, ind)'*(model.d*model.Ainv + model.beta*AinvEETAinv)*model.K_uf(:, ind)...
          -K_fuAinvEMT - K_fuAinvEMT';
    end
    gK_uu = model.d*model.invK_uu ...
            - model.d*model.Ainv/model.beta - AinvEETAinv;
    gBeta = 0;
    gK_ufBase = -(model.d*model.Ainv + model.beta*AinvEETAinv)*model.K_uf ...
        + model.beta*AinvE*M';
    
    for i = 1:length(model.blockEnd)
      ind = gpBlockIndices(model, i);
      invK_uuK_ufDinv = model.invK_uu*model.K_uf(:, ind)*model.Dinv{i};
      gK_uu = gK_uu + model.beta*invK_uuK_ufDinv*blockQ{i}*invK_uuK_ufDinv';
      
      gK_uf(:, ind) = (gK_ufBase(:, ind) ...
                       -model.beta*invK_uuK_ufDinv*blockQ{i})*model.Dinv{i};
      
      
      g_Lambda{i} = 0.5*model.Dinv{i}*blockQ{i}*model.Dinv{i}*model.beta;
      gBeta = gBeta - sum(diag((g_Lambda{i})))/(model.beta*model.beta);
    end
    gK_uu = gK_uu*0.5;
    fhandle = str2func([model.betaTransform 'Transform']);
    gBeta = gBeta*fhandle(model.beta, 'gradfact');
  else
    gK_uu = zeros(model.k, model.k);
    gK_uf = zeros(model.k, model.N);
    for i = 1:length(model.blockEnd)
      if i == 1
        indLen = model.blockEnd(1);
      else
        indLen = model.blockEnd(i) - model.blockEnd(i-1);
      end
      g_Lambda{i} = zeros(indLen, indLen);
    end
    gBeta = 0;
    for j = 1:model.d
      e = zeros(model.k, 1);
      for i = 1:length(model.blockEnd)
        ind = gpDataIndices(model, j, i);
        e = e + model.K_uf(:, ind)*model.Dinv{i, j}*M(ind, j);
      end
      Ainve = model.Ainv{j}*e;
      AinveeT = Ainve*e';
      AinveeTAinv = AinveeT*model.Ainv{j};
      for i = 1:length(model.blockEnd)
        ind = gpDataIndices(model, j, i);
        K_fuAinveyT = model.beta*model.K_uf(:, ind)'*Ainve*M(ind, j)';
        blockQ{i} = -model.D{i, j} + model.beta*M(ind, j)*M(ind, j)' ...
            + model.K_uf(:, ind)'*(model.Ainv{j} + model.beta*AinveeTAinv)*model.K_uf(:, ind)...
            -K_fuAinveyT - K_fuAinveyT';
      end
      gK_uu = gK_uu + model.invK_uu ...
            - model.Ainv{j}/model.beta - AinveeTAinv;
    
      for i = 1:length(model.blockEnd)
        ind = gpDataIndices(model, j, i);
        gK_ufBase = -(model.Ainv{i} + model.beta*AinveeTAinv)*model.K_uf(:, ind) ...
            + model.beta*Ainve*M(ind, j)';
        invK_uuK_ufDinv = model.invK_uu*model.K_uf(:, ind)*model.Dinv{i,j};
        gK_uu = gK_uu + model.beta*invK_uuK_ufDinv*blockQ{i}*invK_uuK_ufDinv';
      
        gK_uf(:, ind) = gK_uf(:, ind) + (gK_ufBase ...
                         -model.beta*invK_uuK_ufDinv*blockQ{i})*model.Dinv{i,j};
      
        if i == 1
          localInd = ind;
        else
          localInd = ind - (model.blockEnd(i-1));
        end
        g_Lambda{i}(localInd, localInd) = g_Lambda{i}(localInd, localInd) ...
            + 0.5*model.Dinv{i,j}*blockQ{i}*model.Dinv{i,j}*model.beta;
      end
    end
    for i = 1:length(model.blockEnd)
      gBeta = gBeta - sum(diag((g_Lambda{i})))/(model.beta*model.beta);
    end
    gK_uu = gK_uu*0.5;
    fhandle = str2func([model.betaTransform 'Transform']);
    gBeta = gBeta*fhandle(model.beta, 'gradfact');
  end
 otherwise
  error('Unknown approximation type');
end
