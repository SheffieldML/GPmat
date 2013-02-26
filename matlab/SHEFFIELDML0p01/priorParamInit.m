function prior = priorParamInit(prior)

% PRIORPARAMINIT Prior model's parameter initialisation.
%
%	Description:
%
%	PRIOR = PRIORPARAMINIT(PRIOR) wrapper function for initialising
%	prior distributions parameters.
%	 Returns:
%	  PRIOR - initialised prior structure.
%	 Arguments:
%	  PRIOR - structure to initialise.
%	
%
%	See also
%	PRIORCREATE


%	Copyright (c) 2003 Neil D. Lawrence


prior = feval([prior.type 'PriorParamInit'], prior);
