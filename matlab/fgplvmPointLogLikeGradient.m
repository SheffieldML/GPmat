function g = fgplvmPointLogLikeGradient(model, x, y)

% FGPLVMPOINTLOGLIKEGRADIENT Log-likelihood gradient for of a point of the GP-LVM.
% FORMAT
% DESC returns the gradient of the log likelihood with respect to
% the latent position, where the log likelihood is conditioned on
% the training set. 
% ARG model : the model for which the gradient computation is being
% done.
% ARG x : the latent position where the gradient is being computed.
% ARG y : the position in data space for which the computation is
% being done.
% RETURN g : the gradient of the log likelihood, conditioned on the
% training data, with respect to the latent position.
%
% SEEALSO : fgplvmPointLogLikelihood, fgplvmOptimisePoint, fgplvmSequenceLogLikeGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM

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
nu(find(isnan(dlnZ_dmu))) = 0.0;
dlnZ_dmu(find(isnan(dlnZ_dmu)))=0.0;

dlnZ_dmu = dlnZ_dmu.*nu;
dlnZ_dvs = 0.5*(dlnZ_dmu.*dlnZ_dmu - nu);

g = dlnZ_dmu*dmu' + dlnZ_dvs*dvs';

if isfield(model, 'prior')
  g = g + priorGradient(model.prior, x);
end
