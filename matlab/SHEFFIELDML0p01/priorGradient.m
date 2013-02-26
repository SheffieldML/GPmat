function g = priorGradient(prior, params)

% PRIORGRADIENT Gradient of the prior with respect to its variables
%
%	Description:
%
%	G = PRIORGRADIENT(PRIOR, PARAMS) wrapper function that computes the
%	gradient of the prior with respect to its parameters.
%	 Returns:
%	  G - gradients of the prior with respect to the parameters.
%	 Arguments:
%	  PRIOR - the prior structure for which the gradients are being
%	   computed.
%	  PARAMS - the parameter locations for which the gradients are being
%	   computed.
%	
%
%	See also
%	PRIORCREATE


%	Copyright (c) 2003 Neil D. Lawrence

  
fhandle = str2func([prior.type 'PriorGradient']);
g = fhandle(prior, params);