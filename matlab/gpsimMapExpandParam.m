function model = gpsimMapExpandParam(model, params)

% GPSIMMAPEXPANDPARAM Expand the given parameters into a GPSIMMAP structure.
% FORMAT
% DESC takes the given vector of parameters and places them in the
% model structure, it then updates any stored representations that
% are dependent on those parameters, for example kernel matrices
% etc..
% ARG model : the model structure for which parameters are to be
% updated.
% ARG params : a vector of parameters for placing in the model
% structure.
% RETURN model : a returned model structure containing the updated
% parameters.
% 
% SEEALSO : gpsimMapCreate, gpsimMapExtractParam, modelExtractParam, gpsimMapUpdateKernels
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFIED : Pei Gao, 2008
  
% SHEFFIELDML
params = real(params);
if isfield(model, 'fix')
  for i = 1:length(model.fix)
    params(model.fix(i).index) = model.fix(i).value;
  end
end

if length(params) ~= model.numParams
  error('Parameter vector is incorrect length');
end
startVal = 1;
endVal = model.kern.nParams;
model.kern = kernExpandParam(model.kern, params(startVal:endVal));

ngParamk = model.ngParam/model.numGenes;

if ~isfield(model, 'isGroupNonlinearity')
  nBaseParam = 3;
else
  if model.isGroupNonlinearity && strcmp(model.nonLinearity{1}, 'activation')
    nBaseParam = 3;
  else 
    nBaseParam = 4;
  end
end

nParamk = nBaseParam + ngParamk;
fhandle = str2func([model.Transform 'Transform']);

for i = 1:model.numGenes
  if isfield(model,'bTransform') && isempty(model.bTransform)
    model.B(i) = params(endVal+nParamk*(i-1)+1);
  else
    model.B(i) = fhandle(params(endVal+nParamk*(i-1)+1), 'atox');
  end
  model.S(i) = fhandle(params(endVal+nParamk*(i-1)+2), 'atox');
  model.D(i) = fhandle(params(endVal+nParamk*(i-1)+3), 'atox');
    
  if nBaseParam == 4
    if isfield(model,'alphaTransform') && isempty(model.alphaTransform)
      model.alpha(i) = params(endVal+nParamk*(i-1)+4);
    else
      model.alpha(i) = fhandle(params(endVal+nParamk*(i-1)+4), 'atox');
    end    
  end   
    
  if model.ngParam > 0
    startV = endVal+nParamk*(i-1)+nBaseParam+1;
    endV = endVal+nParamk*(i-1)+nParamk;
    model.gParam(:,i) =  fhandle(params(startV:endV), 'atox');
  end
end   

if isfield(model, 'includeNoise') && model.includeNoise
  noiseSd = fhandle(params((end-model.numGenes+1):end), 'atox');
  model.noiseVar = noiseSd.*noiseSd;
end

% Check if there is a mean function.
if isfield(model, 'meanFunction') & ~isempty(model.meanFunction)
  startVal = endVal + 1;
  endVal = endVal + model.meanFunction.numParams;
  model.meanFunction = modelExpandParam(model.meanFunction, ...
                                        params(startVal:endVal));
end

model = gpsimMapUpdateKernels(model);
model = gpsimMapUpdatePosteriorCovariance(model);
model.updateW = true; % Needs to be re-done after parameter change.
