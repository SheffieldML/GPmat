function params = gpExtractParam(model)

% GPEXTRACTPARAM Extract a parameter vector from a GP model.
%
% params = gpExtractParam(model)
%

% Copyright (c) 2006 Neil D. Lawrence
% gpExtractParam.m version 1.3



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
 case 'nftc'
  fhandle = str2func([model.betaTransform 'Transform']);
  params =  [kernExtractParam(model.kern) scaleParams fhandle(model.beta, 'xtoa')];
end

