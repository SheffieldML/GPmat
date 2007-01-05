function model = kbrParamInit(model)

% KBRPARAMINIT KBR model parameter initialisation.
% FORMAT
% DESC initialises the kernel based regression
%  model structure with some default parameters.
% ARG model : the model structure which requires initialisation.
% RETURN model : the model structure with the default parameters placed in.
%
% SEEALSO : kbrCreate, modelCreate, modelParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

model.A = randn(model.numData, outputDim)/sqrt(model.numData+1);
model.bias = randn(1, outputDim)/sqrt(model.numData+1);



