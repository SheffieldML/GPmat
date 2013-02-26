function L = orderedNoiseLikelihood(noise, mu, varsigma, y)

% ORDEREDNOISELIKELIHOOD Likelihood of the data under the ORDERED noise model.
%
%	Description:
%
%	ORDEREDNOISELIKELIHOOD(NOISE, MU, VARSIGMA, Y) returns the
%	likelihood of a data set under the  ordered categorical noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the likelihood is required.
%	  MU - input mean locations for the likelihood.
%	  VARSIGMA - input variance locations for the likelihood.
%	  Y - target locations for the likelihood.
%	
%
%	See also
%	ORDEREDNOISEPARAMINIT, ORDEREDNOISELOGLIKELIHOOD, NOISELIKELIHOOD


%	Copyright (c) 2004, 2005 Neil D. Lawrence



D = size(y, 2);

L = zeros(size(mu));
L(find(isnan(y))) = 1;

c = 1./sqrt(noise.variance + varsigma);
for j = 1:D
  % Do lowest category first
  index = find(y(:, j) == 0);
  if ~isempty(index)
    mu(index, j) = mu(index, j) + noise.bias(j) ;
    L(index, j) = cumGaussian(-mu(index, j).*c(index, j));
  end
  
  % Intermediate categories
  index = find(y(:, j) > 0 & y(:, j) < noise.C-1);
  if ~isempty(index)
    for i = index'
      mu(i, j) = mu(i, j) + noise.bias(j) - sum(noise.widths(1:y(i, j)-1));
      L(i, j) = cumGaussian(mu(i, j).*c(i, j)) ...
          - cumGaussian((mu(i, j) ...
                         - noise.widths(y(i, j))) ... 
                        .* c(i, j));
    end
  end
  % Highest category
  index = find(y(:, j) == noise.C-1);
  if ~isempty(index)
    for i = index'
      mu(i, j) = mu(i, j) + noise.bias(j) - sum(noise.widths(1:y(i, j)-1));
    end
    L(index, j) = cumGaussian(mu(index, j).*c(index, j));
  end
end
