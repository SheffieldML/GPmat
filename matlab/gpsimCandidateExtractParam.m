function [param, names] = gpsimCandidateExtractParam(model)

% GPSIMCANDIDATEEXTRACTPARAM Extract the parameters of a GPSIM model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a Gaussian process for single input motif
% modelling.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
%
% SEEALSO : gpsimCreate, gpsimAddCandidate, gpsimCandidateExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2007

% SHEFFIELDML


if nargout>1
  [param, names] = kernExtractParam(model.candidate.kern);
  for i=1:length(model.candidate.mu);
    names{end+1}=['Basal transcription ' num2str(i)];
  end
else
  param = kernExtractParam(model.candidate.kern);
end
fhandle = str2func([model.candidate.bTransform 'Transform']);
param = [param fhandle(model.candidate.B, 'xtoa')];


if isfield(model, 'fix')
  for i = 1:length(model.candidate.fix)
    param(model.candidate.fix(i).index) = model.candidate.fix(i).value;
  end
end
param = real(param);

% Remove main kernel parameters.
for i = model.kern.nParams:-1:1
  if nargout>1
    names(i)=[];
  end
  param(i) = [];
end
