function params = kbrExtractParam(model);

% KBREXTRACTPARAM Extract weights from a kernel based regression model.

params = [model.A(:)' model.b];