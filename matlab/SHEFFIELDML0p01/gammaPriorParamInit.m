function prior = gammaPriorParamInit(prior)

% GAMMAPRIORPARAMINIT Gamma prior model's parameter initialisation.
%
%	Description:
%
%	PRIOR = GAMMAPRIORPARAMINIT(PRIOR) initialises the parameters of the
%	gamma prior with some default parameters.
%	 Returns:
%	  PRIOR - prior structure with initial values in place.
%	 Arguments:
%	  PRIOR - prior structure to be initialised.
%	
%
%	See also
%	PRIORCREATE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


prior.a = 1e-6;
prior.b = 1e-6;

prior.transforms.index = [1 2];
prior.transforms.type = optimiDefaultConstraint('positive');
prior.nParams = 2;
