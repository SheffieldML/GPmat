function H = modelHessian(params, model, varargin)

% MODELHESSIAN Hessian of error function to minimise for given model.
% FORMAT
% DESC gives the Hessian of the objective function for a model. By
% default the objective function is a negative log likelihood.
% ARG params : parameter vector to evaluate at.
% ARG model : model structure to optimise.
% ARG P1, P2, P3 ... : optional additional arguments.
% RETURN H : the Hessian of the error function to be minimised.
% 
% SEEALSO : modelLogLikeHessian, modelObjective, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

fhandle = [model.type 'Hessian'];
if exist(fhandle) == 2
  fhandle = str2func(fhandle);
  H = fhandle(params, model, varargin{:});
else
  fhandle = str2func([model.type 'ExpandParam']);
  model = fhandle(model, params);
  fhandle = str2func([model.type 'LogLikeHessian']);
  H = - fhandle(model, varargin{:});
end
