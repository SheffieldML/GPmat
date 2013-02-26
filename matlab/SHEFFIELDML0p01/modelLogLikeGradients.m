function g = modelLogLikeGradients(model)

% MODELLOGLIKEGRADIENTS Compute a model's gradients wrt log likelihood.
%
%	Description:
%
%	G = MODELLOGLIKEGRADIENTS(MODEL) is a wrapper function to compute
%	the gradients of the log likelihood of a given model.
%	 Returns:
%	  G - teh gradients of the likelihood with respect to the
%	   parameters.
%	 Arguments:
%	  MODEL - the model for which likelihoods are computed.
%	
%
%	See also
%	MODELCREATE


%	Copyright (c) 2006, 2005 Neil D. Lawrence


fhandle = str2func([model.type 'LogLikeGradients']);
g = fhandle(model);

if isfield(model, 'paramGroups')
  g = g*model.paramGroups;
end
