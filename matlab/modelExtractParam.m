function params = modelExtractParam(model)

% MODELEXTRACTPARAM Extract the parameters of a model.

% MLTOOLS

fhandle = str2func([model.type 'ExtractParam']);
params = fhandle(model);