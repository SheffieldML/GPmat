function [mu, varsigma] = dnetPosteriorMeanVar(model, X);

% DNETPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
%
%	Description:
%
%	[MU, SIGMA] = DNETPOSTERIORMEANVAR(MODEL, X) returns the posterior
%	mean and variance for a given set of points.
%	 Returns:
%	  MU - the mean of the posterior distribution.
%	  SIGMA - the variances of the posterior distributions.
%	 Arguments:
%	  MODEL - the model for which the posterior will be computed.
%	  X - the input positions for which the posterior will be computed.
%	
%
%	See also
%	DNETCREATE


%	Copyright (c) 2008 Neil D. Lawrence


mu = dnetOut(model, X);

% Return magnification factors instead of variance.
g = modelOutputGradX(model.mapping, X);
varsigma = zeros(size(X, 1), 1);
for n = 1:size(X, 1)
  gTemp = squeeze(g(n, :, :)); 
  varsigma(n) = exp(0.5*logdet(gTemp*gTemp'));
end

%varsigma = repmat(1/model.beta, size(X, 1), 1);
