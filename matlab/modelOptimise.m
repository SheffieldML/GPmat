function model = modelOptimise(model)

% MODELOPTIMISE Optimise the given model.

% MLTOOLS

fhandle = str2func([model.type 'Optimise']);
model = fhandle(model);