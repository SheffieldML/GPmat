function model = gpsimCandidateExpandParam(model, params)

% GPSIMCANDIDATEEXPANDPARAM Expand the given parameters for a candidate gene.
% FORMAT
% DESC takes the given vector of parameters for a new candidate gene and places them in the
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
% SEEALSO : gpsimExpandParam, gpsimCandidateExtractParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2007

% SHEFFIELDML
  
params = real(params);
if isfield(model, 'fix')
  for i = 1:length(model.candidate.fix)
    params(model.candidate.fix(i).index) = model.candidate.fix(i).value;
  end
end

if length(params) ~= model.candidate.numParams
  error('Parameter vector is incorrect length');
end
baseKernParam = kernExtractParam(model.kern);
startVal = 1;
endVal = model.candidate.kern.nParams-model.kern.nParams;
% pad out parameters for candidate kernel with those from true kernel.
model.candidate.kern = kernExpandParam(model.candidate.kern, [baseKernParam ...
                    params(startVal:endVal)]);

fhandle = str2func([model.candidate.bTransform 'Transform']);
model.candidate.B = fhandle(params(endVal+1:end), 'atox');

% The decays and sensitivities are actually stored in the kernel.
% We'll put them here as well for convenience.
% Only take ones associated with candidate kernel.
counter = 0;
for i = model.kern.numBlocks+1:model.candidate.kern.numBlocks
  counter = counter + 1;
  model.candidate.D(counter) = model.candidate.kern.comp{i}.decay;
  model.candidate.S(counter) = sqrt(model.candidate.kern.comp{i}.variance);
end
model.candidate.mu = model.candidate.B./model.candidate.D;

model = gpsimCandidateUpdateKernels(model);
lengthObs = size(model.candidate.t, 1);
ind = 1:lengthObs;
for i = 1:model.candidate.numGenes
  model.candidate.m(ind) = model.candidate.y(ind) - model.candidate.mu(i)*ones(lengthObs, 1);
  ind = ind + lengthObs;
end
