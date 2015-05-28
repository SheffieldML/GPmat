function f = gpObjective(params, model)

% GPOBJECTIVE Wrapper function for GP objective.
% FORMAT
% DESC returns the negative log likelihood of a Gaussian process
% model given the model structure and a vector of parameters. This
% allows the use of NETLAB minimisation functions to find the model
% parameters.
% ARG params : the parameters of the model for which the objective
% will be evaluated.
% ARG model : the model structure for which the objective will be
% evaluated.
% RETURN f : the negative log likelihood of the GP model.
%
% SEEALSO : scg, conjgrad, gpCreate, gpGradient, gpLogLikelihood, gpOptimise
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% GP

model = gpExpandParam(model, params);
f = - gpLogLikelihood(model);
