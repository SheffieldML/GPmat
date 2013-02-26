function err = modelObjective(params, model, varargin)

% MODELOBJECTIVE Objective function to minimise for given model.
%
%	Description:
%
%	ERR = MODELOBJECTIVE(PARAMS, MODEL, ...) gives the objective
%	function for a model. By default it is the negative log likelihood.
%	 Returns:
%	  ERR - the error function to be minimised.
%	 Arguments:
%	  PARAMS - parameter vector to evaluate at.
%	  MODEL - model structure to optimise.
%	  ... - optional additional arguments.
%	
%
%	See also
%	MODELLOGLIKELIHOOD, MODELGRADIENT, MODELOPTIMISE


%	Copyright (c) 2006 Neil D. Lawrence



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
