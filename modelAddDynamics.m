function model = modelAddDynamics(model, type, varargin)

% MODELADDDYNAMICS Add a dynamics kernel to the model.
% FORMAT
% DESC adds a dynamics model to a model.
% ARG model : the model to add dynamics to.
% ARG type : the type of dynamics model to add in.
% ARG P1, P2, P3, ... : additional arguments to be passed on creation of
% the dynamics model.
%
% SEEALSO : modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2007

% MLTOOLS

type = [type 'Dynamics'];
model.dynamics = modelCreate(type, model.q, model.q, model.X, varargin{:});
params = modelExtractParam(model);
model = modelExpandParam(model, params);

