function params = modelExtractParam(model)

% MODELEXTRACTPARAM Extract the parameters of a model.

fhandle = str2func([model.type 'ExtractParam']);
params = fhandle(model);