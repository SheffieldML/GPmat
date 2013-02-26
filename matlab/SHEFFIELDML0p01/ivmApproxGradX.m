function g = ivmApproxGradX(model, x, m, beta);

% IVMAPPROXGRADX Returns the gradient of the approxmate log-likelihood wrt x.
%
%	Description:
%
%	G = IVMAPPROXGRADX(MODEL, X, M, BETA) returns the gradient of the
%	approximate log likelihood for the IVM with respect to a given input
%	position.
%	 Returns:
%	  G - the gradient of the log likelihood with respect to x.
%	 Arguments:
%	  MODEL - the model for which the gradient is being computed.
%	  X - the input location for which the gradient is to be computed.
%	  M - the output position where the gradient is to be computed.
%	  BETA - the output noise level for which the gradient is being
%	   computed.
%	
%
%	See also
%	IVMPOSTERIORMEANVAR, IVMPOSTERIORGRADMEANVAR


%	Copyright (c) 2004, 2005 Neil D. Lawrence


[mu, varsigma] = ivmPosteriorMeanVar(model, x);
[dmu, dvs] = ivmPosteriorGradMeanVar(model, x);


D = size(m, 2);
nu = 1./(1./beta+varsigma);
dlnZ_dmu = zeros(size(nu));
for i = 1:D
  dlnZ_dmu(:, i) = m(:, i) - mu(:, i) - model.noise.bias(i);
end
dlnZ_dmu = dlnZ_dmu.*nu;
dlnZ_dvs = 0.5*(dlnZ_dmu.*dlnZ_dmu - nu);

g = dlnZ_dmu*dmu' + dlnZ_dvs*dvs';
