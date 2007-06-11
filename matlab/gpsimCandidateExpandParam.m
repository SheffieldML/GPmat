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

% GPSIM
  
model.candidate = gpsimExpandParam(model.candidate, params);
model = gpsimCandidateUpdateKernels(model);