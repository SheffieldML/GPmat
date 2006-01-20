function params = linearExtractParam(model);

% LINEAREXTRACTPARAM Extract weights from a linear model.

% MLTOOLS

params = [model.W(:)' model.b];