function model = gpsimMapUpdateYpredVar (model)
  
% GPSIMMAPUPDATEYPREDVAR Update the variance for y.
% FORMAT
% DESC updates the variance for y. 
% ARG model : the model for which the variance for y is to be
% updated.
% RETURN model : the model with the variance of y.
%
% SEEALSO : gpsimMapCreate, gpsimMapUpdeteYpred
  
% COPYRIGHT : Pei Gao, Magnus Rattray and Neil D. Lawrence, 2008

% SHEFFIELDML
  
  iter = 100;

  for i = 1:iter
    predModel.samp{i} = model;
    f = model.f + sqrt(model.varf).* ...
        randn(size(model.f));
%     figure(1); hold on;
%     plot(model.mapt, f);
    predModel.samp{i} = gpsimMapFunctionalExpandParam(model, f');
%     figure(2); hold on;
%     plot(model.mapt, predModel.samp{i}.ypred(:,1)); 
    
    for j = 1:model.numGenes
      ypred{j}(:,i) = predModel.samp{i}.ypred(:,j);
    end
  end

  for j = 1:model.numGenes
    if isfield(model,'includeNoise') && model.includeNoise
      ypredVar(:,j) = var(ypred{j},0,2) + model.noiseVar(j);
    else
      ypredVar(:,j) = var(ypred{j},0,2);
    end
  end
  
  model.ypredVar = ypredVar;
%  plot(model.mapt, model.ypred(:,1)-2*sqrt(ypredVar(:,1)), 'r-');
