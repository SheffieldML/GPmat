function prior = uniformPriorParamInit(prior)

% UNIFORMPRIORPARAMINIT Uniform prior model's parameter initialisation.
%
%	Description:
%
%	PRIOR = UNIFORMPRIORPARAMINIT(PRIOR) initialises the parameters of
%	the 2 parameter uniform prior with some default parameters.
%	 Returns:
%	  PRIOR - prior structure with initial values in place.
%	 Arguments:
%	  PRIOR - prior structure to be initialised.
%	
%
%	See also
%	PRIORCREATE


%	Copyright (c) 2012 Antti Honkela


prior.a = 0;
prior.b = 1;

prior.nParams = 0;
prior.isBounded = 1;
