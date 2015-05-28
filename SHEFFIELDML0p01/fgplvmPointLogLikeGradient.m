function g = fgplvmPointLogLikeGradient(model, x, y)

% FGPLVMPOINTLOGLIKEGRADIENT Log-likelihood gradient for of a point of the GP-LVM.
%
%	Description:
%
%	G = FGPLVMPOINTLOGLIKEGRADIENT(MODEL, X, Y) returns the gradient of
%	the log likelihood with respect to the latent position, where the
%	log likelihood is conditioned on the training set.
%	 Returns:
%	  G - the gradient of the log likelihood, conditioned on the
%	   training data, with respect to the latent position.
%	 Arguments:
%	  MODEL - the model for which the gradient computation is being
%	   done.
%	  X - the latent position where the gradient is being computed.
%	  Y - the position in data space for which the computation is being
%	   done.
%	
%
%	See also
%	FGPLVMPOINTLOGLIKELIHOOD, FGPLVMOPTIMISEPOINT, FGPLVMSEQUENCELOGLIKEGRADIENT


%	Copyright (c) 2006 Neil D. Lawrence


logTwoPi = log(2*pi);
[mu, varSigma] = gpPosteriorMeanVar(model, x);
[dmu, dvs] = gpPosteriorGradMeanVar(model, x);


% For more general models this should be done with the noise
% toolbox (see ivmGradX in the ivm toolbox for more details).
nu = 1./varSigma;
dlnZ_dmu = zeros(size(nu));
for i = 1:model.d
  dlnZ_dmu(:, i) = y(:, i) - mu(:, i);
end
dlnZ_dmu = dlnZ_dmu.*nu;
dlnZ_dvs = 0.5*(dlnZ_dmu.*dlnZ_dmu - nu);

g = dlnZ_dmu*dmu' + dlnZ_dvs*dvs';

if isfield(model, 'prior')
  g = g + priorGradient(model.prior, x);
end
