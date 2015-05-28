function model = multimodelParamInit(model)

% MULTIMODELPARAMINIT MULTIMODEL model parameter initialisation.
%
%	Description:
%
%	MODEL = MULTIMODELPARAMINIT(MODEL) initialises the multi-task
%	learning wrapper model structure with some default parameters.
%	 Returns:
%	  MODEL - the model structure with the default parameters placed in.
%	 Arguments:
%	  MODEL - the model structure which requires initialisation.
%	
%
%	See also
%	MULTIMODELCREATE, MODELCREATE, MODELPARAMINIT


%	Copyright (c) 2007 Neil D. Lawrence


for i = 1:length(model.comp)
  model.comp{i} = modelParamInit(model.comp{i});
end