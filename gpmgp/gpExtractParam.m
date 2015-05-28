function [params, names] = gpExtractParam(model)

% GPEXTRACTPARAM Extract a parameter vector from a GP model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a Gaussian process.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
%
% DESC does the same as above, but also returns parameter names.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
% RETURN names : cell array of parameter names.
%
% SEEALSO : gpCreate, gpExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2007, 2009

% GP

if nargout > 1
  returnNames = true;
else
  returnNames = false;
end

% Check if the output scales are being learnt.
if model.learnScales
  fhandle = str2func([model.scaleTransform 'Transform']);
  scaleParams = fhandle(model.scale, 'xtoa');
  if returnNames
    for i = 1:length(scaleParams)
      scaleParamNames{i} = ['Output Scale ' num2str(i)];
    end
  end
else
  scaleParams = [];
  scaleParamNames = {};
end

% Check if there is a mean function.
if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
  if returnNames
    [meanFuncParams, meanFuncParamNames] = modelExtractParam(model.meanFunction);
    for i = 1:length(meanFuncParamNames)
      meanFuncParamNames{i} = ['Mean Func, ' meanFuncParamNames{i}];
    end
  else
    meanFuncParams = modelExtractParam(model.meanFunction);
  end
else
  meanFuncParamNames = {};
  meanFuncParams =[];
end
if returnNames
  [kernParams, kernParamNames] = kernExtractParam(model.kern);
  for i = 1:length(kernParamNames)
    kernParamNames{i} = ['Kernel, ' kernParamNames{i}];
  end
else
  kernParams = kernExtractParam(model.kern);
end
switch model.approx
 case 'ftc'
  params =  [kernParams meanFuncParams scaleParams];
  if returnNames
    names = {kernParamNames{:}, meanFuncParamNames{:}, scaleParamNames{:}};
  end
  if model.optimiseBeta
    fhandle = str2func([model.betaTransform 'Transform']);
    betaParam = fhandle(model.beta, 'xtoa');
    params = [params betaParam(:)'];
    if returnNames
      for i = 1:length(betaParam)
        betaParamNames{i} = ['Beta ' num2str(i)];
      end
      names = {names{:}, betaParamNames{:}};
    end
  end
 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
  paramPart = [kernParams meanFuncParams scaleParams];
  if returnNames
    names = {kernParamNames{:}, meanFuncParamNames{:}, ...
             scaleParamNames{:}};
  end
  if model.optimiseBeta
    fhandle = str2func([model.betaTransform 'Transform']);
    betaParam = fhandle(model.beta, 'xtoa');
    paramPart = [paramPart betaParam(:)'];
    if returnNames
      for i = 1:length(betaParam)
        betaParamNames{i} = ['Beta ' num2str(i)];
      end
      names = {names{:}, betaParamNames{:}};
    end
  end
  if model.fixInducing
    params = paramPart;
  else
    params =  [model.X_u(:)' paramPart];
    if returnNames 
      for i = 1:size(model.X_u, 1)
        for j = 1:size(model.X_u, 2)
          X_uNames{i, j} = ['X_u(' num2str(i) ', ' num2str(j) ')'];
        end
      end
      names = {X_uNames{:}, names{:}};
    end
  end
end
