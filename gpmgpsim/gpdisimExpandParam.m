function model = gpdisimExpandParam(model, params)

% GPDISIMEXPANDPARAM Expand the given parameters into a GPDISIM structure.
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
if model.includeNoise,
  simMultiKern = model.kern.comp{1};
else
  simMultiKern = model.kern;
end

model.delta = simMultiKern.comp{2}.di_decay;
model.sigma = sqrt(simMultiKern.comp{2}.di_variance);
for i = 2:simMultiKern.numBlocks
  model.D(i-1) = simMultiKern.comp{i}.decay;
  model.S(i-1) = sqrt(simMultiKern.comp{i}.variance);
end
model.mu = model.B./model.D;

model = gpsimUpdateKernels(model);
lengthObs = size(model.t, 1);
ind = 1:lengthObs;
model.m(ind) = model.y(ind);
ind = ind + lengthObs;
for i = 1:model.numGenes
  model.m(ind) = model.y(ind) - model.mu(i)*ones(lengthObs, 1);
  ind = ind + lengthObs;
end
