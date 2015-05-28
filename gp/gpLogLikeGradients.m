function [gParam, gX_u, gX, g_beta] = gpLogLikeGradients(model, X, M, X_u)

% GPLOGLIKEGRADIENTS Compute the gradients for the parameters and X.
% FORMAT
% DESC computes the gradients of the Gaussian process log
% likelihood with respect to the parameters of the model.
% ARG : model : the model structure for which gradients are computed.
% RETURN gParam : the gradient of the log likelihood with respect to
% the model parameters.
%
% DESC computes the gradients of the Gaussian process log
% likelihood with respect to the parameters of the model and with
% respect to any inducing variables.
% ARG : model : the model structure for which gradients are computed.
% RETURN gParam : the gradient of the log likelihood with respect to
% the model parameters.
% RETURN gX_u : the gradient of the log likelihood with respect to
% the inducing variables. If inducing variables aren't being used
% this returns zero.
%
% DESC computes the gradients of the Gaussian process log
% likelihood with respect to the parameters of the model, with
% respect to any inducing variables and with respect to input data
% locations. This is used for computing gradients in the GP-LVM.
% ARG : model : the model structure for which gradients are computed.
% RETURN gParam : the gradient of the log likelihood with respect to
% the model parameters (including any gradients with respect to beta).
% RETURN gX_u : the gradient of the log likelihood with respect to
% the inducing variables. If inducing variables aren't being used
% this returns zero.
% RETURN gX : the gradient of the log likelihood with respect to
% the input data locations.
%
% DESC computes the gradients of the Gaussian process log
% likelihood with respect to the parameters of the model, with
% respect to any inducing variables and with respect to input data
% locations. This is used for computing gradients in the GP-LVM.
% ARG : model : the model structure for which gradients are computed.
% RETURN gParam : the gradient of the log likelihood with respect to
% the model parameters.
% RETURN gX_u : the gradient of the log likelihood with respect to
% the inducing variables. If inducing variables aren't being used
% this returns zero.
% RETURN gX : the gradient of the log likelihood with respect to
% the input data locations.
% RETURN gbeta : the gradient of the log likelihood with respect to beta.
%
% DESC computes the gradients of the Gaussian process log
% likelihood with respect to the model parameters (and optionally,
% as above with respect to inducing variables and input data) given
% the target data, input data and inducing variable
% locations. 
% ARG : model : the model structure for which gradients are computed.
% ARG : X : the input data locations for which gradients are computed.
% ARG : M : the scaled and bias removed target data for which the 
% gradients are computed.
% ARG : X_U : the inducing variable locations for which gradients are computed.
% RETURN gParam : the gradient of the log likelihood with respect to
% the model parameters.
%
% SEEALSO : gpLogLikelihood, modelLogLikeGradients, fgplvmLogLikeGradients
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2007, 2009
%
% MODIFICATIONS : Carl Henrik Ek, 2008

% GP

  if nargin < 4
    if isfield(model, 'X_u')
      X_u = model.X_u;
    else
      X_u = [];
    end
    if nargin < 3 && ~isfield(model, 'S')
      M = model.m;
    end
    if nargin < 2
      X = model.X;
    end
  end

  gX_u = [];
  gX = [];

  g_scaleBias = gpScaleBiasGradient(model);
  if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
    g_meanFunc = gpMeanFunctionGradient(model);
  else
    g_meanFunc = [];
  end

  switch model.approx
   case 'ftc'
    % Full training conditional.
    if nargout > 2
      %%% Prepare to Compute Gradients with respect to X %%%
      gKX = kernGradX(model.kern, X, X);
      gKX = gKX*2;
      dgKX = kernDiagGradX(model.kern, X);
      for i = 1:model.N
        gKX(i, :, i) = dgKX(i, :);
      end
      gX = zeros(model.N, model.q);
    end
    
    %%% Gradients of Kernel Parameters %%%
    g_param = zeros(1, model.kern.nParams);
    if isfield(model, 'beta')
      g_beta = 0;
    else
      g_beta = [];
    end

    % For very high D, we use the matrix S which is M*M'
    if isfield(model, 'S')
      gK = localSCovarianceGradients(model);
      if nargout > 2
        %%% Compute Gradients with respect to X %%%
        counter = 0;
        for i = 1:model.N
          counter = counter + 1;
          for j = 1:model.q
            gX(i, j) = gX(i, j) + gKX(:, j, i)'*gK(:, counter);
          end
        end
      end
      %%% Compute Gradients of Kernel Parameters %%%
      g_param = g_param + kernGradient(model.kern, X, gK);
    else
      for k = 1:model.d
        gK = localCovarianceGradients(model, M(:, k), k);
        if nargout > 2
          %%% Compute Gradients with respect to X %%%
          ind = gpDataIndices(model, k);
          counter = 0;
          for i = ind
            counter = counter + 1;
            for j = 1:model.q
              gX(i, j) = gX(i, j) + gKX(ind, j, i)'*gK(:, counter);
            end
          end
        end
        %%% Compute Gradients of Kernel Parameters %%%
        if model.isMissingData
          g_param = g_param ...
                    + kernGradient(model.kern, ...
                                   X(model.indexPresent{k}, :), ...
                                   gK);
        else
            g_param = g_param + kernGradient(model.kern, X, gK);
        end
      end
      if isfield(model, 'beta') && model.optimiseBeta
        if size(model.beta, 1) == 1
          g_beta = g_beta + sum(diag(gK));
        elseif size(model.beta, 2)==1 ...
              & size(model.beta, 1)==model.N
          g_beta = g_beta + diag(gK);
        elseif size(model.beta, 2) == model.d ...
              & size(model.beta, 1) == model.N
          g_beta(:, k) = diag(gK);
        else
          error('Unusual dimensions for model.beta.');
        end
      end
    end
    
   case {'dtc', 'dtcvar', 'fitc', 'pitc'}
    % Sparse approximations.
    [gK_u, gK_uf, gK_star, g_beta] = gpCovGrads(model, M);
    
    %%% Compute Gradients of Kernel Parameters %%%
    gParam_u = kernGradient(model.kern, X_u, gK_u);
    gParam_uf = kernGradient(model.kern, X_u, X, gK_uf);

    g_param = gParam_u + gParam_uf;
    
    %%% Compute Gradients with respect to X_u %%%
    gKX = kernGradX(model.kern, X_u, X_u);
    
    % The 2 accounts for the fact that covGrad is symmetric
    gKX = gKX*2;
    dgKX = kernDiagGradX(model.kern, X_u);
    for i = 1:model.k
      gKX(i, :, i) = dgKX(i, :);
    end
    
    if ~model.fixInducing | nargout > 1
      % Allocate space for gX_u
      gX_u = zeros(model.k, model.q);
      % Compute portion associated with gK_u
      for i = 1:model.k
        for j = 1:model.q
          gX_u(i, j) = gKX(:, j, i)'*gK_u(:, i);
        end
      end
      
      % Compute portion associated with gK_uf
      gKX_uf = kernGradX(model.kern, X_u, X);
      for i = 1:model.k
        for j = 1:model.q
          gX_u(i, j) = gX_u(i, j) + gKX_uf(:, j, i)'*gK_uf(i, :)';
        end
      end
    end
    if nargout > 2
      %%% Compute gradients with respect to X %%%
      
      % Allocate space for gX
      gX = zeros(model.N, model.q);
      
      % this needs to be recomputed so that it is wrt X not X_u
      gKX_uf = kernGradX(model.kern, X, X_u);
      
      for i = 1:model.N
        for j = 1:model.q
          gX(i, j) = gKX_uf(:, j, i)'*gK_uf(:, i);
        end
      end    
    end
   otherwise
    error('Unknown model approximation.')
  end


  switch model.approx
   case 'ftc'
    % Full training conditional. Nothing required here.
   case 'dtc'
    % Deterministic training conditional.  
   case {'fitc', 'dtcvar'}
    % Fully independent training conditional.
    % Variational sparse approximation.
    
    if nargout > 2
      % deal with diagonal term's effect on X gradients..
      gKXdiag = kernDiagGradX(model.kern, X);
      for i = 1:model.N
        gX(i, :) = gX(i, :) + gKXdiag(i, :)*gK_star(i);
      end
    end
    
    % deal with diagonal term's affect on kernel parameters.
    g_param = g_param + kernDiagGradient(model.kern, X, gK_star);


   case 'pitc'
    % Partially independent training conditional.
    
    if nargout > 2
      % deal with block diagonal term's effect on X gradients.
      startVal = 1;
      for i = 1:length(model.blockEnd)
        endVal = model.blockEnd(i);
        ind = startVal:endVal;
        gKXblock = kernGradX(model.kern, X(ind, :), X(ind, :));
        
        % The 2 accounts for the fact that covGrad is symmetric
        gKXblock = gKXblock*2;
        
        % fix diagonal
        dgKXblock = kernDiagGradX(model.kern, X(ind, :));
        for j = 1:length(ind)
          gKXblock(j, :, j) = dgKXblock(j, :);
        end
        
        for j = ind
          for k = 1:model.q
            subInd = j - startVal + 1;
            gX(j, k) = gX(j, k) + gKXblock(:, k, subInd)'*gK_star{i}(:, subInd);
          end
        end
        startVal = endVal + 1;
      end
    end
    % deal with block diagonal's effect on kernel parameters.
    for i = 1:length(model.blockEnd);
      ind = gpBlockIndices(model, i);
      g_param = g_param ...
                + kernGradient(model.kern, X(ind, :), gK_star{i});
    end
    
   otherwise
    error('Unrecognised model approximation');
  end

  if nargout < 4
    if (~isfield(model, 'optimiseBeta') && ~strcmp(model.approx, 'ftc')) ...
          | model.optimiseBeta
      % append beta gradient to end of parameters
      gParam = [g_param(:)' g_meanFunc g_scaleBias g_beta];
    else
      gParam = [g_param(:)' g_meanFunc g_scaleBias];
    end
  else
    gParam = [g_param(:)' g_meanFunc g_scaleBias];
  end

  % if there is only one output argument, pack gX_u and gParam into it.
  if nargout == 1;
    gParam = [gX_u(:)' gParam];
  end
end

function gK = localCovarianceGradients(model, y, dimension)

% LOCALCOVARIANCEGRADIENTS

  if ~isfield(model, 'isSpherical') || model.isSpherical
    invKy = model.invK_uu*y;
    gK = -model.invK_uu + invKy*invKy';
  else
    if model.isMissingData
      m = y(model.indexPresent{dimension});
    else
      m = y;
    end
    invKy = model.invK_uu{dimension}*m;
    gK = -model.invK_uu{dimension} + invKy*invKy';
  end
  gK = gK*.5;
end    

function gK = localSCovarianceGradients(model)

% LOCALCOVARIANCEGRADIENTS

  gK = -model.d*model.invK_uu + model.invK_uu*model.S*model.invK_uu;
  gK = gK*.5;
end  
