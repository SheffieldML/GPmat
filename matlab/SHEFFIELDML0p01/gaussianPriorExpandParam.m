function prior = gaussianPriorExpandParam(prior, params)

% GAUSSIANPRIOREXPANDPARAM Expand Gaussian prior structure from param vector.
%
%	Description:
%
%	PRIOR = GAUSSIANPRIOREXPANDPARAM(PRIOR, PARAMS) places the given
%	parameters in a Gaussian prior structure.
%	 Returns:
%	  PRIOR - the structure with the parameters in place.
%	 Arguments:
%	  PRIOR - the structure to place the parameters in.
%	  PARAMS - the parameters to place in the structure.
%	
%
%	See also
%	% SEEALSO PRIOREXPANDPARAM


%	Copyright (c) 2004 Neil D. Lawrence


prior.precision = params(1);
