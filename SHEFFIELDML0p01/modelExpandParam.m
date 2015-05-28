function model = modelExpandParam(model, params, dim)

% MODELEXPANDPARAM Update a model structure with parameters.
%
%	Description:
%
%	MODEL = MODELEXPANDPARAM(MODEL, PARAM) returns a model structure
%	filled with the parameters in the given vector. This is used as a
%	helper function to enable parameters to be optimised in, for
%	example, the NETLAB optimisation functions.
%	 Returns:
%	  MODEL - model structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  MODEL - the model structure in which the parameters are to be
%	   placed.
%	  PARAM - vector of parameters which are to be placed in the model
%	   structure.
%	
%	
%
%	See also
%	MODELEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2005, 2006 Neil D. Lawrence


%	With modifications by Cark Henrik Ek 2007


if isfield(model, 'paramGroups')
  params = params*model.paramGroups';
end

fhandle = str2func([model.type 'ExpandParam']);
if(nargin<3)
  model = fhandle(model, params);
else
  model = fhandle(model,params,dim);
end