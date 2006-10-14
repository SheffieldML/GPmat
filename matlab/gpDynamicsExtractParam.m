function param = gpDynamicsExtractParam(model)

% GPDYNAMICSEXTRACTPARAM Extract parameters from the GP dynamics model.

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
