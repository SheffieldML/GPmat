function x = optimiMinimize(objectiveGradient, params, options, varargin);

% OPTIMIMINIMIZE Wrapper for Carl Rasmussen's minimize function.
% FORMAT
% DESC is a nother wrapper for Carl Rasmussen's minimize function,
% but this time using functions which return the gradient and the
% function value together (such as gpObjectiveGradient). 
% ARG objectiveGradient : function that returns the value (as a
% ARG params : the initial parameters to be optimised.
% scalar) and the gradients (as a row vector) of the function to be
% optimised.
% ARG options : a NETLAB style options vector.
% ARG P1, P2, P3 ... : further arguments to be passed to
% OBJECTIVEGRADIENT.
% RETURN params : row vector of optimised parameters.
%
% SEEALSO : minimize, cgcarl
%
% COPYRIGHT : Neil D. Lawrence, 2006

% OPTIMI

x = minimize(params', objectiveGradient, options(14), ...
             varargin{:})';

