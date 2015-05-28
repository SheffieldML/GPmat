function [f, g] = gradFuncWrapper(x, func, grad, varargin)

% GRADFUNCWRAPPER Wrapper function to enable use of Carl Rasmussen's minimze function.

% OPTIMI

func = str2func(func);
grad = str2func(grad);
f = func(x', varargin{:});

if nargout > 1
  g = grad(x', varargin{:})';
end
