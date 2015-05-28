function model = gpsimExpandParam(model, params)

% GPSIMEXPANDPARAM Expand the given parameters into a GPSIM structure.
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
% SEEALSO : gpsimCreate, gpsimExtractParam, modelExtractParam, gpsimUpdateKernels
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

fhandle = str2func([model.bTransform 'Transform']);
model.B = fhandle(params(endVal+1:end), 'atox');

% The decays and sensitivities are actually stored in the kernel.
% We'll put them here as well for convenience.
if isfield(model, 'proteinPrior') && ~isempty(model.proteinPrior)
  for i = 2:model.kern.numBlocks
    if model.includeNoise
      model.D(i-1) = model.kern.comp{1}.comp{i}.decay;
      model.S(i-1) = sqrt(model.kern.comp{1}.comp{i}.variance);
    else
      model.D(i-1) = model.kern.comp{i}.decay;
      model.S(i-1) = sqrt(model.kern.comp{i}.variance);
    end
  end  
else
  for i = 1:model.kern.numBlocks
    if model.includeNoise
      model.D(i) = model.kern.comp{1}.comp{i}.decay;
      model.S(i) = sqrt(model.kern.comp{1}.comp{i}.variance);
    else
      model.D(i) = model.kern.comp{i}.decay;
      model.S(i) = sqrt(model.kern.comp{i}.variance);
    end
  end
  
end
model.mu = model.B./model.D;

model = gpsimUpdateKernels(model);

if isfield(model, 'proteinPrior') && ~isempty(model.proteinPrior)
  if model.includeNoise
    yInd = 1:model.kern.comp{1}.diagBlockDim{2};
    mInd = model.kern.comp{1}.diagBlockDim{1} + yInd;
    for i = 1:model.numGenes
      model.m(mInd) = model.y(yInd) - model.mu(i)* ...
          ones(model.kern.comp{1}.diagBlockDim{i+1},1);
      yInd = yInd + model.kern.comp{1}.diagBlockDim{i+1};
      mInd = mInd + model.kern.comp{1}.diagBlockDim{i+1};
    end
  else
    yInd = 1:model.kern.diagBlockDim{2};
    mInd = model.kern.diagBlockDim{1} + yInd;
    for i = 1:model.numGenes
      model.m(mInd) = model.y(yInd) - model.mu(i)*ones(model.kern.diagBlockDim{i+1}, ...
                                                     1);
      yInd = yInd + model.kern.diagBlockDim{i+1};
      mInd = mInd + model.kern.diagBlockDim{i+1};
    end
  end
else
  lengthObs = size(model.t, 1);
  ind = 1:lengthObs;
  for i = 1:model.numGenes
    model.m(ind) = model.y(ind) - model.mu(i)*ones(lengthObs, 1);
    ind = ind + lengthObs;
  end
end
