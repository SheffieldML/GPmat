function f = gpdisimObjective(params, model)

% GPDISIMOBJECTIVE Wrapper function for GPDISIM objective.
% FORMAT
% DESC returns the negative log likelihood of a Gaussian process
% model for single input motifs given the model structure and
% a vector parameters. This allows the use of NETLAB minimisation
% functions to find the model parameters.
% ARG params : the parameters of the model for which the objective
% will be evaluated.
% ARG model : the model structure for which the objective will be
% evaluated.
% RETURN f : the negative log likelihood of the GPDISIM model.
%
% SEEALSO : scg, conjgrad, gpdisimCreate, gpdisimGradient, gpdisimLogLikelihood, gpdisimOptimise
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MODIFICATIONS : Antti Honkela, 2007

% SHEFFIELDML

model = gpdisimExpandParam(model, params);
f = - gpdisimLogLikelihood(model);
