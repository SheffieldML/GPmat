function [x, options, flog, pointlog] = cgcarl(f, x, options, gradf, varargin)

% CGCARL Wrapper for Carl Rasmussen's conjugate gradient implemntation.
%
%	Description:
%	[x, options, flog, pointlog] = cgcarl(f, x, options, gradf, varargin)
%


length = options(14);
x = minimize(x', 'gradFuncWrapper', length, f, gradf, varargin{:})';










