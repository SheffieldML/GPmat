function model = modelExpandParam(model, params)

% MODELEXPANDPARAM Update a model structure with parameters.

fhandle = str2func([model.type 'ExpandParam']);
model = fhandle(model, params);