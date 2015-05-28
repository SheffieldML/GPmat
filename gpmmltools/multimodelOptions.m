function options = multimodelOptions(modelType, number, varargin);

% MULTIMODELOPTIONS Create a default options structure for the MULTIMODEL model.
% FORMAT
% DESC creates a default options structure for the multi-task learning wrapper model
% structure.
% ARG modelType : the model type that the multi-task model is based on.
% ARG number : the number of components in the multi-task model.
% ARG p1, p2, p3, ... : optional additional arguments to be passed to the
% 'sub-model's options structure.
% RETURN options : the default options structure.
%
% SEEALSO : multimodelCreate, modelOptions
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

options.type = modelType;
options.numModels = number;
fhandle = str2func([modelType 'Options']);
options.compOptions = fhandle(varargin{:});
