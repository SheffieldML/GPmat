function g = lfmLogLikeGradients(model)

% LFMLOGLIKEGRADIENTS Compute the gradients of the log likelihood of a LFM model.
%
%	Description:
%
%	G = LFMLOGLIKEGRADIENTS(MODEL) computes the gradients of the log
%	likelihood of the given Gaussian process for use in a single input
%	motif protein network.
%	 Returns:
%	  G - the gradients of the parameters of the model.
%	 Arguments:
%	  MODEL - the model for which the log likelihood is computed.
%	
%
%	See also
%	LFMCREATE, LFMLOGLIKELIHOOD, LFMGRADIENT


%	Copyright (c) 2006 Neil D. Lawrence


covGrad = -model.invK + model.invK*model.m*model.m'*model.invK;
covGrad = 0.5*covGrad;
g = kernGradient(model.kern, model.t, covGrad);



if isfield(model, 'fix')
  for i = 1:length(model.fix)
    g(model.fix(i).index) = 0;
  end
end
