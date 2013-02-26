function [param, names] = lfmExtractParam(model)

% LFMEXTRACTPARAM Extract the parameters of an LFM model.
%
%	Description:
%
%	PARAMS = LFMEXTRACTPARAM(MODEL) extracts the model parameters from a
%	structure containing the information about a Gaussian process latent
%	force model.
%	 Returns:
%	  PARAMS - a vector of parameters from the model.
%	 Arguments:
%	  MODEL - the model structure containing the information about the
%	   model.
%	
%
%	See also
%	LFMCREATE, LFMEXPANDPARAM, MODELEXTRACTPARAM


%	Copyright (c) 2007 Neil D. Lawrence


if nargout>1
  [param, names] = kernExtractParam(model.kern);
else
  param = kernExtractParam(model.kern);
end

if isfield(model, 'fix')
  for i = 1:length(model.fix)
    param(model.fix(i).index) = model.fix(i).value;
  end
end
param = real(param);