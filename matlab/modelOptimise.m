function model = modelOptimise(model)

% MODELOPTIMISE Optimise the given model.

fhandle = str2func([model.type 'Optimise']);
model = fhandle(model);