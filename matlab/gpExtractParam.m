function params = gpExtractParam(model)

% GPEXTRACTPARAM Extract a parameter vector from a GP model.

% FGPLVM

if model.learnScales
  fhandle = str2func([model.scaleTransform 'Transform']);
  scaleParams = fhandle(model.scale, 'xtoa');
else
  scaleParams = [];
end
switch model.approx
 case 'ftc'
  params =  [kernExtractParam(model.kern) scaleParams];
 case {'dtc', 'fitc', 'pitc'}
  fhandle = str2func([model.betaTransform 'Transform']);
  params =  [model.X_u(:)' kernExtractParam(model.kern) ...
             scaleParams fhandle(model.beta, 'xtoa')];
end

