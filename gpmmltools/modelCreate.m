function model = modelCreate(type, numIn, numOut, varargin)

% MODELCREATE Create a model of the specified type.
% FORMAT
% DESC creates a model of the given type.
% ARG type : the type of the model to create, for example, 'kbr'
% for kernel based regression, 'mlp' for multi-layer perceptron,
% 'linear' for a linear model.
% ARG numIn : number of inputs to the model (or latent dimensions
% for latent variable models.
% ARG numOut : number of outputs from the model (or data dimensions
% for latent variable models. 
% ARG P1, P2, P3,... : optional arguments to be passed to the model
% creation code.
% RETURN model : the model created.
%
% SEEALSO : modelExpandParam, modelExtractParam
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% MLTOOLS

fhandle = str2func([type, 'Create']);
model = fhandle(numIn, numOut, varargin{:});
