function g = gpMeanFunctionGradient(model)

% GPMEANFUNCTIONGRADIENT Compute the log likelihood gradient wrt the scales.
% FORMAT
% DESC computes the gradient of the log likelihood with respect to the
% scales. In the future the gradients with respect to the biases
% may also be included. 
% ARG model : the model for which the gradients are to be computed.
% ARG g : the gradients of the likelihood with respect to the
% mean function's parameters.
%
% SEEALSO : gpCreate, gpScaleBiasGradient, gpLogLikeGradients, gpLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2009

% GP

if isfield(model, 'isSpherical') & ~model.isSpherical
  error('Currently only implemented for spherical');
else
  if model.isMissingData
    error('Currently not implemented for missing data.');
  end
end
  
if isfield(model, 'meanFunction') & ~isempty(model.meanFunction)
  g = zeros(1, model.meanFunction.numParams);
  % compute gradients here.
  switch model.approx
   case 'ftc'
    gmu = model.invK_uu*model.m;
   case {'dtc', 'dtcvar'}
    gmu = (model.m - model.K_uf'*model.Ainv*(model.K_uf*model.m))*model.beta;
   case 'fitc'
    Dinvm = model.Dinv*model.m;
    gmu = (Dinvm-(model.Dinv*model.K_uf')...
           *(model.Ainv*model.K_uf)*Dinvm)*model.beta;   
   case 'pitc'
    % Loop through the blocks computing each part to be added.
    gmu = zeros(model.N, model.d);
    K_ufDinvm = zeros(model.k, model.d);
    K_ufDinv = zeros(model.k, model.N);
    for i = 1:length(model.blockEnd)
      ind = gpBlockIndices(model, i);
      Dinvm{i} = model.Dinv{i}*model.m(ind, :);
      K_ufDinvm = K_ufDinvm + model.K_uf(:, ind)*Dinvm{i};
    end
    for i = 1:length(model.blockEnd)
      ind = gpBlockIndices(model, i);
      gmu(ind, :) = (Dinvm{i} ...
                     - model.Dinv{i}...
                     *model.K_uf(:, ind)'*(model.Ainv*K_ufDinvm))*model.beta;
    end
  end
  gmu = gmu./repmat(model.scale, model.N, 1);
  goutputDparam = modelOutputGrad(model.meanFunction, model.X);
  for i = 1:model.meanFunction.numParams
    g(1, i) = sum(sum(gmu.*squeeze(goutputDparam(:, i, :))));
  end
else
  g = [];
end
