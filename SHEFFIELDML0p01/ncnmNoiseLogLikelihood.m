function L = ncnmNoiseLogLikelihood(noise, mu, varsigma, y)

% NCNMNOISELOGLIKELIHOOD Log-likelihood of data under null category noise model.
%
%	Description:
%
%	L = NCNMNOISELOGLIKELIHOOD(NOISE, MU, SIGMA, Y) returns the
%	log-likelihood of a set of targets under a given NCNM noise model
%	with process mean and variances provided.
%	 Returns:
%	  L - scalar value containing the log-likelihood of the targets.
%	 Arguments:
%	  NOISE - the noise model structure.
%	  MU - the mean input into the noise model.
%	  SIGMA - the variance input into the noise model.
%	  Y - the targets whose likelihoods are being predicted.
%	
%
%	See also
%	NOISELOGLIKELIHOOD, NCNMLIKELIHOOD


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


D = size(y, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end

L = 0;
c = 1./sqrt(noise.sigma2 + varsigma);
for j = 1:D
  % Do negative class first.
  index = find(y(:, j) == -1);
  if ~isempty(index)
    mu(index, j) = mu(index, j) + noise.width/2;
    L = L + sum(lnCumGaussian(-mu(index, j).*c(index, j)));
    % account for observations of missingness variable.
    L = L + log((1-noise.gamman))*length(index);
  end
  
  % The null category.
  index = find(isnan(y(:, j)));
  if ~isempty(index)
    mu(index, j) = mu(index, j) + noise.width/2;
    u = mu(index, j).*c(index, j);
    uprime = (mu(index, j) - noise.width).* c(index, j); 
    L = L + sum(lnCumGaussSum(-u, uprime, noise.gamman, noise.gammap));
  end
  % The positive class.
  index = find(y(:, j) == 1);
  if ~isempty(index)
    mu(index, j) = mu(index, j) - noise.width/2;
    L = L + sum(lnCumGaussian(mu(index, j).*c(index, j)));
    % account for observations of missingness variable.
    L = L + log((1-noise.gammap))*length(index);
  end
end
  
