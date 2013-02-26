function prior = gaussianPriorParamInit(prior)

% GAUSSIANPRIORPARAMINIT Gaussian prior model's parameter initialisation.
%
%	Description:
%
%	PRIOR = GAUSSIANPRIORPARAMINIT(PRIOR) initialises the parameters of
%	the Gaussian prior with some default parameters.
%	 Returns:
%	  PRIOR - prior structure with initial values in place.
%	 Arguments:
%	  PRIOR - prior structure to be initialised.
%	
%
%	See also
%	PRIORCREATE, GAMMAPRIORPARAMINIT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


prior.precision = 1;

prior.transforms.index = [1];
prior.transforms.type = optimiDefaultConstraint('positive');
prior.nParams = 1;
