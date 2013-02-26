function options = modelOptions(modelType, varargin)

% MODELOPTIONS Returns a default options structure for the given model.
%
%	Description:
%
%	OPTIONS = MODELOPTIONS(MODELTYPE, ...) returns a default options
%	structure for the given model.
%	 Returns:
%	  OPTIONS - options structure.
%	 Arguments:
%	  MODELTYPE - the type of model.
%	  ... - optional additional arguments (dependent on model type).
%	
%
%	See also
%	MODELCREATE


%	Copyright (c) 2006 Neil D. Lawrence


fhandle = str2func([modelType 'Options']);
options = fhandle(varargin{:});