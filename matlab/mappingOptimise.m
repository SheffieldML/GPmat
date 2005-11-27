function model = mappingOptimise(model, X, Y, varargin)

% MAPPINGOPTIMISE Optimise the given model.

fhandle = str2func([model.type 'Optimise']);
model = fhandle(model, X, Y, varargin{:});