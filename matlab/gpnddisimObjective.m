function f = gpnddisimObjective(params, model)

% GPNDDISIMOBJECTIVE Wrapper function for GPNDDISIM objective.
% FORMAT
% DESC returns the negative log likelihood of a Gaussian process
% model for single input motifs given the model structure and
% a vector parameters. This allows the use of NETLAB minimisation
% functions to find the model parameters.
% ARG params : the parameters of the model for which the objective
% will be evaluated.
% ARG model : the model structure for which the objective will be
% evaluated.
% RETURN f : the negative log likelihood of the GPNDDISIM model.
%
% SEEALSO : scg, conjgrad, gpnddisimCreate, gpnddisimGradient, gpnddisimLogLikelihood, gpnddisimOptimise
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MODIFICATIONS : Antti Honkela, 2007
%
% COPYRIGHT : Jaakko Peltonen, 2011

% GPSIM

try,
    model = gpnddisimExpandParam(model, params);
    f = - gpnddisimLogLikelihood(model);
catch
    f = -Inf;
end
