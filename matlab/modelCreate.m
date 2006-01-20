function model = modelCreate(type, numIn, numOut, varargin)

% MODELCREATE Create a model of the specified type.

% MLTOOLS

fhandle = str2func([type, 'Create']);
model = fhandle(numIn, numOut, varargin{:});