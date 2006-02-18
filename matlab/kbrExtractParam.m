function params = kbrExtractParam(model);

% KBREXTRACTPARAM Extract weights from a kernel based regression model.

% MLTOOLS

params = [model.A(:)' model.bias];