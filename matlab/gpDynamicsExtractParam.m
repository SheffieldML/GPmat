function param = gpDynamicsExtractParam(model)

% GPDYNAMICSEXTRACTPARAM Extract parameters from the GP dynamics model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a Gaussian process dynamics model.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
%
% SEEALSO : gpExtractParam, gpDynamicsCreate, gpDynamicsExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM

param = gpExtractParam(model);

if ~model.learn 
  % If we aren't learning model parameters extract only X_u;
  if ~model.learnScales
    param = param(1:model.k*model.q);
  else
    switch model.approx
     case 'ftc'
      param =  [param(end-model.d + 1:end)];
     case {'dtc', 'fitc', 'pitc'}
      param =  [param(1:model.k*model.q) param(end-model.d:end-1)];
    end
  end
end
