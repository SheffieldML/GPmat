function g = modelGradient(params, model, varargin)

% MODELGRADIENT Gradient of error function to minimise for given model.
% FORMAT
% DESC gives the gradient of the objective function for a model. By
% default the objective function is a negative log likelihood.
% ARG params : parameter vector to evaluate at.
% ARG model : model structure to optimise.
% ARG P1, P2, P3 ... : optional additional arguments.
% RETURN g : the gradient of the error function to be minimised.
% 
% SEEALSO : modelLogLikeGradient, modelObjective, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

fhandle = [model.type 'Gradient'];
if exist(fhandle) == 2
  fhandle = str2func(fhandle);
  g = fhandle(params, model, varargin{:});
else
  fhandle = str2func([model.type 'ExpandParam']);
  model = fhandle(model, params);
  fhandle = str2func([model.type 'LogLikeGradients']);
  g = - fhandle(model, varargin{:});
end
