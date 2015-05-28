function options = modelOptions(modelType, varargin)

% MODELOPTIONS Returns a default options structure for the given model.
% FORMAT
% DESC returns a default options structure for the given model.
% ARG modelType : the type of model.
% ARG P1, P2, P3 ... : optional additional arguments (dependent on
% model type).
% RETURN options : options structure.
%
% SEEALSO : modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

fhandle = str2func([modelType 'Options']);
options = fhandle(varargin{:});
