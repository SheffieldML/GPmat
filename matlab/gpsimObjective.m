function f = gpsimObjective(params, model)

% GPSIMOBJECTIVE Wrapper function for GPSIM objective.
% FORMAT
% DESC returns the negative log likelihood of a Gaussian process
% model for single input motifs given the model structure and
% a vector parameters. This allows the use of NETLAB minimisation
% functions to find the model parameters.
% ARG params : the parameters of the model for which the objective
% will be evaluated.
% ARG model : the model structure for which the objective will be
% evaluated.
% RETURN f : the negative log likelihood of the GPSIM model.
%
% SEEALSO : scg, conjgrad, gpsimCreate, gpsimGradient, gpsimLogLikelihood, gpsimOptimise
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% SHEFFIELDML

model = gpsimExpandParam(model, params);
f = - gpsimLogLikelihood(model);
