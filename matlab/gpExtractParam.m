function params = gpExtractParam(model)

% GPEXTRACTPARAM Extract a parameter vector from a GP model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a Gaussian process.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
%
% SEEALSO : gpCreate, gpExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

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
  paramPart = [kernExtractParam(model.kern) ...
               scaleParams fhandle(model.beta, 'xtoa')];
  if model.fixInducing
    params = paramPart;
  else
    params =  [model.X_u(:)' paramPart];
  end
end

