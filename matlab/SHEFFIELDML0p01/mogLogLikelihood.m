function ll = mogLogLikelihood(model)

% MOGLOGLIKELIHOOD Mixture of Gaussian's log likelihood.
%
%	Description:
%
%	LLL = MOGLOGLIKELIHOOD(MODEL) computes the variational lower bound
%	on the log likelihood of a mixtures of probabilistic PCA model, it
%	wraps the mogLowerBound command.
%	 Returns:
%	  LLL - the lower bound on the log likelihood computed for the
%	   model.
%	 Arguments:
%	  MODEL - the model for which log likelihood is to be computed.
%	
%
%	See also
%	MOGCREATE, MOGLOWERBOUND


%	Copyright (c) 2006, 2008 Neil D. Lawrence

model = mogEstep(model);
ll = mogLowerBound(model);

