function f = dnetObjective(params, model)

% DNETOBJECTIVE Wrapper function for Density Network objective.
% FORMAT
% DESC provides a wrapper function for the Density Network, it
% takes the negative of the log likelihood, feeding the parameters
% correctly to the model.
% ARG params : the parameters of the Density Network model.
% ARG model : the model structure in which the parameters are to be
% placed.
% RETURN f : the negative of the log likelihood of the model.
% 
% SEEALSO : dnetCreate, dnetLogLikelihood, dnetExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

model = dnetExpandParam(model, params);
f = - dnetLogLikelihood(model);
