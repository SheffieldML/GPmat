function model = fgplvmAddDynamics(model, type, varargin)

% FGPLVMADDDYNAMICS Add a dynamics kernel to the model.
%
% model = fgplvmAddDynamics(model, type, varargin)
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmAddDynamics.m version 1.2



type = [type 'Dynamics'];
model.dynamics = modelCreate(type, model.q, model.q, model.X, varargin{:});
params = fgplvmExtractParam(model);
model = fgplvmExpandParam(model, params);

