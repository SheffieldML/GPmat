function model = modelParamInit(model, varargin)

% MODELPARAMINIT Initialise the parameters of the model.
% FORMAT
% DESC initialises the parameters of the model with some sensible
% values.
% ARG model : model for which initialisation will be performed.
% ARG P1, P2, P3, ... : optional additional arguments.
% RETURN model : model with parameters initialised.
%
% SEEALSO : modelCreate
% 
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

fhandle = str2func([model.type 'ParamInit']);
options = fhandle(model, varargin{:});

