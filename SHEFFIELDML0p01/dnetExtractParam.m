function [params, names] = dnetExtractParam(model)

% DNETEXTRACTPARAM Extract weights and biases from an DNET.
%
%	Description:
%
%	[PARAMS, NAMES] = DNETEXTRACTPARAM(MODEL) returns a vector of all
%	the parameters from a density network.
%	 Returns:
%	  PARAMS - vector of all the weights and biases returned by the
%	   model. The structure is governed by dnetpak.
%	  NAMES - optional additional returned cell array of the names of
%	   the parameters.
%	 Arguments:
%	  MODEL - the model from which we wish to extract the weights and
%	   biases.
%	
%
%	See also
%	DNETCREATE, DNETEXPANDPARAM, MODELEXTRACTPARAM


%	Copyright (c) 2008 Neil D. Lawrence


if nargout > 1
  [params, names] = modelExtractParam(model.mapping);
  names{end+1} = 'Output Precision';
else
  params = modelExtractParam(model.mapping);
end
func = str2func([model.betaTransform 'Transform']);
params(end+1) = func(model.beta, 'xtoa');
