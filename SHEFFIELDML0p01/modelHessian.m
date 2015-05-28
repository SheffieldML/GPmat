function H = modelHessian(params, model, varargin)

% MODELHESSIAN Hessian of error function to minimise for given model.
%
%	Description:
%
%	H = MODELHESSIAN(PARAMS, MODEL, ...) gives the Hessian of the
%	objective function for a model. By default the objective function is
%	a negative log likelihood.
%	 Returns:
%	  H - the Hessian of the error function to be minimised.
%	 Arguments:
%	  PARAMS - parameter vector to evaluate at.
%	  MODEL - model structure to optimise.
%	  ... - optional additional arguments.
%	
%
%	See also
%	MODELLOGLIKEHESSIAN, MODELOBJECTIVE, MODELOPTIMISE


%	Copyright (c) 2006 Neil D. Lawrence


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
