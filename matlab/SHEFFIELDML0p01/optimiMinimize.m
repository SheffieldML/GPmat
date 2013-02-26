function x = optimiMinimize(objectiveGradient, params, options, varargin);

% OPTIMIMINIMIZE Wrapper for Carl Rasmussen's minimize function.
%
%	Description:
%
%	PARAMS = OPTIMIMINIMIZE(OBJECTIVEGRADIENT, PARAMS, OPTIONS, ...) is
%	a nother wrapper for Carl Rasmussen's minimize function, but this
%	time using functions which return the gradient and the function
%	value together (such as gpObjectiveGradient).
%	 Returns:
%	  PARAMS - row vector of optimised parameters.
%	 Arguments:
%	  OBJECTIVEGRADIENT - function that returns the value (as a
%	  PARAMS - the initial parameters to be optimised. scalar) and the
%	   gradients (as a row vector) of the function to be optimised.
%	  OPTIONS - a NETLAB style options vector.
%	  ... - further arguments to be passed to OBJECTIVEGRADIENT.
%	
%
%	See also
%	MINIMIZE, CGCARL


%	Copyright (c) 2006 Neil D. Lawrence


x = minimize(params', objectiveGradient, options(14), ...
             varargin{:})';

