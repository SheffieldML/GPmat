function model = kbrExpandParam(model, params,dim);

% KBREXPANDPARAM Create model structure from KBR model's parameters.
%
%	Description:
%
%	MODEL = KBREXPANDPARAM(MODEL, PARAMS) takes a vector of KBR weights
%	and centres and places them in their respective positions in the KBR
%	model.
%	 Returns:
%	  MODEL - the model with the weights distributed in the correct
%	   places.
%	 Arguments:
%	  MODEL - the model in which the weights are to be placed.
%	  PARAMS - a vector of the weights to be placed in the model.
%	
%	
%	
%
%	See also
%	KBRUNPAK, KBRCREATE, KBREXTRACTPARAM


%	Copyright (c) 2008 Neil D. Lawrence


%	With modifications by Carl Henrik Ek 2008

if(nargin<3)
  startVal = 1;
  endVal = model.numData*model.outputDim;
  model.A = reshape(params(1:endVal), model.numData, model.outputDim);
  model.bias = params(endVal+1:end);
else
  model.A(:,dim) = params(1:1:model.numData*length(dim));
  model.bias(dim) = params(model.numData*length(dim)+1:1:end);
end