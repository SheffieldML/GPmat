function g = gpSequenceLogLikeGradient(model, X, Y)

% GPSEQUENCELOGLIKEGRADIENT Log-likelihood gradient for of a sequence of the GP-LVM.
% FORMAT
% DESC is a helper function for fgplvmSequenceLogLikeGradient. It returns the gradient of the log likelihood with respect to
% the input positions, where the log likelihood is conditioned on
% the training set. It ignores any priors over the input positions.
% ARG model : the model for which the gradient computation is being
% done.
% ARG X : the latent sequence where the gradient is being computed.
% ARG Y : the sequence in data space for which the computation is
% being done.
% RETURN g : the gradient of the log likelihood, conditioned on the
% training data, with respect to the latent sequence.
%
% SEEALSO : fgplvmSequenceLogLikelihood, fgplvmOptimiseSequence, fgplvmSequenceLogLikeGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM

logTwoPi = log(2*pi);
if model.isMissingData
  [mu, covarSigma] = gpPosteriorMeanCovar(model, X);
  [dmu, dcovar] = gpPosteriorGradMeanCovar(model, X);
else
  [mu, covarSigma, factors] = gpPosteriorMeanCovar(model, X);
  [dmu, dcovar, factors] = gpPosteriorGradMeanCovar(model, X);
end

% For more general models this should be done with the noise
% toolbox (see ivmGradX in the ivm toolbox for more details).
missing = true;
if ~any(isnan(Y))
  U = jitChol(covarSigma);
  missing = false;
end
Ydiff = (Y-mu);
dlnZ_dmu = zeros(size(mu));
for j = 1:model.d
  dlnZ_dcovar{j} = zeros(size(X, 1));
  if missing
    ind = find(~isnan(Ydiff(:, j)));
    if length(ind) ~= 0
      U = jitChol(covarSigma(ind, ind));
    else
      U = [];
    end
  else
    ind = [1:size(Ydiff, 1)]';
  end
  UinvYdiff = U'\Ydiff(ind, j);
  dlnZ_dmu(ind, j) = U\UinvYdiff/factors(j);
  dlnZ_dcovar{j}(ind, ind) = 0.5*(UinvYdiff*UinvYdiff'/factors(j) - eye(length(ind)));
  dlnZ_dcovar{j}(ind, ind) = U\dlnZ_dcovar{j}(ind, ind)/U';
  dlnZ_dcovar{j}(ind, ind) = dlnZ_dcovar{j}(ind, ind)/factors(j);
end

g = zeros(size(X, 1), model.q);
for k = 1:model.d
  for j = 1:model.q
    g(:, j) = g(:, j) + dlnZ_dmu(:, k).*dmu{j}(:, k);
  end
  for j = 1:model.q
    % Since dlnZ_dcovar is symmetric we double dcovar apart from
    % the diagonal term ...
    diagCovar = diag(dcovar{j});
    dcovar2 = 2*dcovar{j};
    dcovar2(1:size(dcovar2, 2)+1:end) = diagCovar;
    for  i = 1:size(X, 1)
      g(i, j) = g(i, j) + dcovar2(:, i)'*dlnZ_dcovar{k}(:, i)*factors(k);
    end
  end
end
  
