function model = fgplvmAddDynamics(model, type, varargin)

% FGPLVMADDDYNAMICS Add a dynamics kernel to the model.
% FORMAT
% DESC adds a dynamics model to the FGPLVM.
% ARG model : the model to add dynamics to.
% ARG type : the type of dynamics model to add in.
% ARG P1, P2, P3, ... : additional arguments to be passed on creation of
% the dynamics model.
%
% SEEALSO : modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2007
%
% MODIFICATIONS : Carl Henrik Ek, 2008
%
% FGPLVM

type = [type 'Dynamics'];
model.dynamics = modelCreate(type, model.q, model.q, model.X, varargin{:});
params = fgplvmExtractParam(model);
model.dynamics.numParams = length(modelExtractParam(model.dynamics));
model.numParams = model.numParams + model.dynamics.numParams;
model = fgplvmExpandParam(model, params);

