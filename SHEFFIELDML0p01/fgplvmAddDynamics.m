function model = fgplvmAddDynamics(model, type, varargin)

% FGPLVMADDDYNAMICS Add a dynamics kernel to the model.
%
%	Description:
%
%	FGPLVMADDDYNAMICS(MODEL, TYPE, ...) adds a dynamics model to the
%	FGPLVM.
%	 Arguments:
%	  MODEL - the model to add dynamics to.
%	  TYPE - the type of dynamics model to add in.
%	  ... - additional arguments to be passed on creation of the
%	   dynamics model.
%	
%	
%	
%
%	See also
%	MODELCREATE


%	Copyright (c) 2005, 2006, 2007 Neil D. Lawrence


%	With modifications by Carl Henrik Ek 2008

type = [type 'Dynamics'];
model.dynamics = modelCreate(type, model.q, model.q, model.X, varargin{:});
params = fgplvmExtractParam(model);
model.dynamics.numParams = length(modelExtractParam(model.dynamics));
model.numParams = model.numParams + model.dynamics.numParams;
model = fgplvmExpandParam(model, params);

