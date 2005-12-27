function params = gpExtractParam(model)

% GPEXTRACTPARAM Extract a parameter vector from a GP model.

% FGPLVM

switch model.approx
 case 'ftc'
  params =  kernExtractParam(model.kern);
 case {'dtc', 'fitc', 'pitc'}
  fhandle = str2func([model.sigma2Transform 'Transform']);
  params =  [model.X_u(:)' kernExtractParam(model.kern) ...
             fhandle(model.sigma2, 'xtoa')];
end

