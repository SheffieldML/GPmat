function params = linearExtractParam(model);

% LINEAREXTRACTPARAM Extract weights from a linear model.

params = [model.W(:)' model.b];