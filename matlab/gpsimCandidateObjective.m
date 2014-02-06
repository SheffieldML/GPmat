function f = gpsimCandidateObjective(params, model)

% GPSIMCANDIDATEOBJECTIVE Wrapper function for GPSIM candidate gene objective.
% FORMAT
% DESC returns the negative log likelihood of candidate genes given a
% Gaussian process  model for single input motifs given the model
% structure and a vector parameters. This allows the use of NETLAB
% minimisation functions to find the model parameters.
% ARG params : the parameters of the model for which the objective
% will be evaluated.
% ARG model : the model structure for which the objective will be
% evaluated.
% RETURN f : the negative log likelihood of the candidate genes given the
% GPSIM model.
%
% SEEALSO : scg, conjgrad, gpsimCreate, gpsimAddCandidate, gpsimCandidateGradient, gpsimCandidateLogLikelihood, gpsimCandidateOptimise
% 
% COPYRIGHT : Neil D. Lawrence, 2007

% SHEFFIELDML

model = gpsimCandidateExpandParam(model, params);
f = - gpsimCandidateLogLikelihood(model);
