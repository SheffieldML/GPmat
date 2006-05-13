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

% GPSIM

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
for i = 1:model.kern.numBlocks
  model.D(i) = model.kern.comp{i}.decay;
  model.S(i) = sqrt(model.kern.comp{i}.variance);
end
model.mu = model.B./model.D;

model = gpsimUpdateKernels(model);
lengthObs = size(model.t, 1);
ind = 1:lengthObs;
for i = 1:model.numGenes
  model.m(ind) = model.y(ind) - model.mu(i)*ones(lengthObs, 1);
  ind = ind + lengthObs;
end
