function params = kbrExtractParam(model,dim);

% KBREXTRACTPARAM Extract parameters from the KBR model structure.
%
%	Description:
%
%	PARAM = KBREXTRACTPARAM(MODEL) extracts parameters from the kernel
%	based regression model structure into a vector of parameters for
%	optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the model.
%	 Arguments:
%	  MODEL - the model structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = KBREXTRACTPARAM(MODEL) extracts parameters and
%	parameter names from the kernel based regression model structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the model.
%	  NAMES - cell array of strings containing names for each parameter.
%	 Arguments:
%	  MODEL - the model structure containing the parameters to be
%	   extracted.
%	
%	
%
%	See also
%	% SEEALSO KBRCREATE, KBREXPANDPARAM, MODELEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2005, 2006 Neil D. Lawrence
% 	kbrExtractParam.m version 1.4


if(nargin<2)
  params = [model.A(:)' model.bias];
else
  params = model.A(:,dim);
  params = [params(:)' model.bias(dim)];
end
