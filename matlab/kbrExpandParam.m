function model = kbrExpandParam(model, params,dim);

% KBREXPANDPARAM Create model structure from KBR model's parameters.
%
%	Description:
%
%	MODEL = KBREXPANDPARAM(MODEL, PARAM) returns a kernel based
%	regression model structure filled with the parameters in the given
%	vector. This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
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
%	See also
%	KBRCREATE, KBREXTRACTPARAM, MODELEXPANDPARAM


%	Copyright (c) 2005, 2006 Neil D. Lawrence
% 	kbrExpandParam.m version 1.4


if(nargin<3)
  startVal = 1;
  endVal = model.numData*model.outputDim;
  model.A = reshape(params(1:endVal), model.numData, model.outputDim);
  model.bias = params(endVal+1:end);
else
  model.A(:,dim) = params(1:1:model.numData*length(dim));
  model.bias(dim) = params(model.numData*length(dim)+1:1:end);
end