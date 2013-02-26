function ll = rbfperiodicLogLikelihood(model)

% RBFPERIODICLOGLIKELIHOOD Log likelihood of RBFPERIODIC model.
%
%	Description:
%
%	LL = RBFPERIODICLOGLIKELIHOOD(MODEL) computes the log likelihood of 
%	the periodic radial basis function model.
%	 Returns:
%	  LL - the computed log likelihood.
%	 Arguments:
%	  MODEL - the model structure for which log likelihood is being
%	   computed.
%	
%
%	See also
%	RBFPERIODICCREATE, RBFPERIODICLOGLIKEGRADIENTS, MODELLOGLIKELIHOOD


%	Copyright (c) 2007 Neil D. Lawrence


ypred = rbfperiodicOut(model, model.X);

ll = - 0.5.*(sum(sum((model.y-ypred).^2)));
ll = real(ll);