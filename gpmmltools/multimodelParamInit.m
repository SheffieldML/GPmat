function model = multimodelParamInit(model)

% MULTIMODELPARAMINIT MULTIMODEL model parameter initialisation.
% FORMAT
% DESC initialises the multi-task learning wrapper
%  model structure with some default parameters.
% ARG model : the model structure which requires initialisation.
% RETURN model : the model structure with the default parameters placed in.
%
% SEEALSO : multimodelCreate, modelCreate, modelParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

for i = 1:length(model.comp)
  model.comp{i} = modelParamInit(model.comp{i});
end
