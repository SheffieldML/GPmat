function model = mappingOptimise(model, X, Y, varargin)

% MAPPINGOPTIMISE Optimise the given model.
%
%	Description:
%	model = mappingOptimise(model, X, Y, varargin)
%

fhandle = str2func([model.type 'Optimise']);
model = fhandle(model, X, Y, varargin{:});