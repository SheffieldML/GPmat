function L = ncnmNoiseLikelihood(noise, mu, varsigma, y)

% NCNMNOISELIKELIHOOD Likelihood of data under null category noise model.
%
%	Description:
%
%	L = NCNMNOISELIKELIHOOD(NOISE, MU, SIGMA, Y) returns the likelihood
%	of a set of targets under a given NCNM noise model with process mean
%	and variances provided.
%	 Returns:
%	  L - vector of a likelihood for each of the targets.
%	 Arguments:
%	  NOISE - the noise model structure.
%	  MU - the mean input into the noise model.
%	  SIGMA - the variance input into the noise model.
%	  Y - the targets whose likelihoods are being predicted.
%	
%
%	See also
%	NOISELIKELIHOOD, NCNMLOGLIKELIHOOD


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


D = size(y, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end


c = 1./sqrt(noise.sigma2 + varsigma);
for j = 1:D
  % Negatively labelled data.
  index = find(y(:, j) == -1);
  if ~isempty(index)
    mu(index, j) = mu(index, j) + noise.width/2;
    L(index, j) = cumGaussian(-mu(index, j).*c(index, j));
    % account for observations of missingness variable.
    L(index,j) = L(index,j)*(1-noise.gamman);
  end
  
  % Missing data.
  index = find(isnan(y(:, j));
  if ~isempty(index)
    mu(index, j) = mu(index, j) + noise.width/2;
    L(index, j) = noise.gamman*cumGaussian(-mu(index, j).*c(index, j)) ...
	+ noise.gammap*cumGaussian((mu(index, j) - noise.width) .* c(index, j));
  end
  % Highest category
  index = find(y(:, j) == 1);
  if ~isempty(index)
    mu(index, j) = mu(index, j) - noise.width/2;
    L(index, j) = cumGaussian(mu(index, j).*c(index, j));
    L(index,j) = L(index,j)*(1-noise.gammap);
  end
end
