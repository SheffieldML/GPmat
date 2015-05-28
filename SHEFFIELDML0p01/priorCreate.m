function prior = priorCreate(type)

% PRIORCREATE Create a prior structure given a type.
%
%	Description:
%
%	PRIOR = PRIORCREATE(TYPE) creates a prior structure given a type.
%	 Returns:
%	  PRIOR - The prior structure.
%	 Arguments:
%	  TYPE - Type of prior to be created,  some standard types are
%	   'gamma', 'gaussian', 'laplace' and 'invgamma'.
%	
%
%	See also
%	PRIORPARAMINIT


%	Copyright (c) 2003 Neil D. Lawrence


prior.type = type;
prior = priorParamInit(prior);
