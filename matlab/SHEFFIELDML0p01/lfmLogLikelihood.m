function ll = lfmLogLikelihood(model)

% LFMLOGLIKELIHOOD Compute the log likelihood of a LFM model.
%
%	Description:
%
%	LL = LFMLOGLIKELIHOOD(MODEL) computes the log likelihood of the
%	given Gaussian process for use in a latent force model.
%	 Returns:
%	  LL - the log likelihood of the data set.
%	 Arguments:
%	  MODEL - the model for which the log likelihood is computed.
%	
%
%	See also
%	LFMCREATE, LFMLOGLIKEGRADIENT, LFMOBJECTIVE


%	Copyright (c) 2007 Neil D. Lawrence


dim = size(model.y, 1);
ll = -dim*log(2*pi) - model.logDetK - model.m'*model.invK*model.m;
ll = ll*0.5;
