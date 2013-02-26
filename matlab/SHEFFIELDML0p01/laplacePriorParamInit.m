function prior = laplacePriorParamInit(prior)

% LAPLACEPRIORPARAMINIT Laplace prior model's parameter initialisation.
%
%	Description:
%
%	PRIOR = LAPLACEPRIORPARAMINIT(PRIOR) initialises the parameters of
%	the Laplace prior with some default parameters.
%	 Returns:
%	  PRIOR - prior structure with initial values in place.
%	 Arguments:
%	  PRIOR - prior structure to be initialised.
%	
%
%	See also
%	PRIORCREATE, GAMMAPRIORPARAMINIT, GAUSSIANPRIORPARAMINIT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


prior.precision = 1;

prior.transforms.index = [1];
prior.transforms.type = optimiDefaultConstraint('positive');
prior.nParams = 1;
