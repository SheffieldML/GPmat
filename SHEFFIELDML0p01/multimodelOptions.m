function options = multimodelOptions(modelType, number, varargin);

% MULTIMODELOPTIONS Create a default options structure for the MULTIMODEL model.
%
%	Description:
%
%	OPTIONS = MULTIMODELOPTIONS(MODELTYPE, NUMBER, ...) creates a
%	default options structure for the multi-task learning wrapper model
%	structure.
%	 Returns:
%	  OPTIONS - the default options structure.
%	 Arguments:
%	  MODELTYPE - the model type that the multi-task model is based on.
%	  NUMBER - the number of components in the multi-task model.
%	  ... - optional additional arguments to be passed to the
%	   'sub-model's options structure.
%	
%
%	See also
%	MULTIMODELCREATE, MODELOPTIONS


%	Copyright (c) 2007 Neil D. Lawrence


options.type = modelType;
options.numModels = number;
fhandle = str2func([modelType 'Options']);
options.compOptions = fhandle(varargin{:});