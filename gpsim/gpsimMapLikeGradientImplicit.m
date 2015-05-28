function [g dlogPdf] = gpsimMapLikeGradientImplicit(model)

% GPSIMMAPLIKEGRADIENTIMPLICIT computes the implicit part of the gradient

% FORMAT
% ARG model : the model for which the data component of the Hessian is
% to be updated.
% RETURN g :
% RETURN dlogPdf : 
%
% SEEALSO : gpsimMapCreate, gpsimMapFunctionalGradient, gpsimMapUpdatePosteriorCovariance
%
% COPYRIGHT : Pei Gao, Magnus Rattray and Neil D. Lawrence, 2008

% SHEFFIELDML
  
g1 = zeros(1, model.kern.nParams);
g2 = zeros(1, model.numParams-model.kern.nParams);

% covGrad = inv(eye(size(model.K))+model.K*model.W); 
covGrad = model.covf*model.invK;
funcGrad = model.invK*model.f;

for ftime_index = 1:model.numMapPts
    kernGrad(ftime_index, :) = kernGradient(model.kern, model.mapt, covGrad)*funcGrad(ftime_index);
%    dfuncGrad(ftime_index,:) = zeros(size(g2));
    dfuncGrad(ftime_index, :) = gpsimMapFunctionalLikeGrad2(model, ftime_index);
end

modelParamGrad = covGrad*model.K*dfuncGrad;

for ftime_index = 1:model.numMapPts
    dWdf(:,:,ftime_index) = gpsimMapFunctionalWGradient(model, ftime_index);
    
    dlogPdf(ftime_index) = -0.5*trace(covGrad*model.K*dWdf(:,:,ftime_index));
    g1 = g1 + dlogPdf(ftime_index)*kernGrad(ftime_index, :);
    
    g2 = g2 + dlogPdf(ftime_index)*modelParamGrad(ftime_index,:);
end

% Check if model parameters are being optimised in a transformed space
[paramvec names] = gpsimMapExtractParam(model);

if isfield(model,'bTransform') && isempty(model.bTransform)
  Bindex = [];
  for i = 1:length(names)
    if strcmp(names{i}(1:5), 'Basal')
      Bindex = [Bindex i];
    end
  end
  paramvec(Bindex) = 0;
end

if isfield(model,'alphaTransform') && isempty(model.alphaTransform)
  alphaIndex = [];
  for i = 1:length(names)
    if strcmp(names{i}(1:5), 'Param')
      alphaIndex = [alphaIndex i];
    end
  end
  paramvec(alphaIndex) = 0;
end

modelParams = paramvec((model.kern.nParams+1):end);

if ~isempty(model.Transform)
  fhandle = str2func([model.Transform 'Transform']);
  modelParams = fhandle(modelParams, 'atox');
  g2 = g2.*fhandle(modelParams, 'gradfact');
end

g = [g1 g2];

