function [param, names] = gpsimMapExtractParam(model)

% GPSIMMAPEXTRACTPARAM Extract the parameters of a GPSIMMAP model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a Gaussian process for single input motif
% modelling.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
%
% SEEALSO : gpsimMapCreate, gpsimMapExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006
%  
% MODIFIED : Pei Gao, 2008

% SHEFFIELDML

if nargout>1
  [param, names] = kernExtractParam(model.kern);
else
  param = kernExtractParam(model.kern);
end

fhandle = str2func([model.Transform 'Transform']);

for i = 1:model.numGenes
  if isfield(model,'bTransform') && isempty(model.bTransform);
    param = [param model.B(i)];
  else
    param = [param fhandle(model.B(i), 'xtoa')];
  end
    
  param = [param fhandle(model.S(i), 'xtoa') fhandle(model.D(i), ...
                                                    'xtoa')];
  
  if  isfield(model, 'includeRepression') && model.includeRepression
    if isfield(model,'alphaTransform') && isempty(model.alphaTransform);
      param = [param model.alpha(i)];
    else
      param = [param fhandle(model.alpha(i), 'xtoa')];
    end   
  end
    
  if model.ngParam > 0
    param = [param (fhandle(model.gParam(:,i), 'xtoa'))'];   
  end
    
  if nargout >1
    names{end+1} = ['Basal' num2str(i)];
    names{end+1} = ['Sensitivity' num2str(i)];
    names{end+1} = ['Decay' num2str(i)];
    
    if isfield(model, 'isGroupNonlinearity') && model.isGroupNonlinearity
      if strcmp(model.nonLinearity{i}, 'repression')
        names{end+1} = ['alpha' num2str(i)];
      end
    else      
      names{end+1} = ['alpha' num2str(i)]; 
    end        
         
    if model.ngParam > 0
      ngParamk = model.ngParam/model.numGenes;
      names{(end+1):(end+ngParamk)} = ['Parameters for g' num2str(i)];
    end
  end
end

if isfield(model,'includeNoise') && model.includeNoise
  param = [param fhandle(sqrt(model.noiseVar), 'xtoa')];
  if nargout > 1
    for i=1:model.numGenes
      names = {names{:} ['noiseSd' num2str(i)]};
    end
  end
end

if isfield(model, 'fix')
  for i = 1:length(model.fix)
    param(model.fix(i).index) = model.fix(i).value;
  end
end

param = real(param);

% Check if there is a mean function.
if isfield(model, 'meanFunction') & ~isempty(model.meanFunction)
  if nargout>1
    [meanFuncParams, meanFuncNames] = modelExtractParam(model.meanFunction);
  else
    meanFuncParams = modelExtractParam(model.meanFunction);
  end
else
  meanFuncParams =[];
  meanFuncNames = {};
end

param = [param meanFuncParams];
if nargout > 1
  names = {names{:} meanFuncNames{:}};
end
