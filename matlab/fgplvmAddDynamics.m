function model = fgplvmAddDynamics(model, type, varargin)

% FGPLVMADDDYNAMICS Add a dynamics kernel to the model.

% FGPLVM

type = [type 'Dynamics'];
model.dynamics = modelCreate(type, model.q, model.q, model.X, varargin{:});
params = fgplvmExtractParam(model);
model = fgplvmExpandParam(model, params);

