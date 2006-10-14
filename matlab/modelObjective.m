function err = modelObjective(params, model, varargin)

% MODELOBJECTIVE Objective function to minimise for given model.
% FORMAT
% DESC gives the objective function for a model. By default it is
% the negative log likelihood.
% ARG params : parameter vector to evaluate at.
% ARG model : model structure to optimise.
% ARG P1, P2, P3 ... : optional additional arguments.
% RETURN err : the error function to be minimised.
% 
% SEEALSO : modelLogLikelihood, modelGradient, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS


fhandle = [model.type 'Objective'];
if exist(fhandle) == 2
  fhandle = str2func(fhandle);
  err = fhandle(params, model, varargin{:});
else
  fhandle = str2func([model.type 'ExpandParam']);
  model = fhandle(model, params);
  fhandle = str2func([model.type 'LogLikelihood']);
  err = - fhandle(model, varargin{:});
end
