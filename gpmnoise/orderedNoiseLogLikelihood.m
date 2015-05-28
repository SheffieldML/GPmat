function L = orderedNoiseLogLikelihood(noise, mu, varsigma, y)


% ORDEREDNOISELOGLIKELIHOOD Log likelihood of the data under the ORDERED noise model.
% FORMAT
% DESC returns the log likelihood of a data set under the  ordered categorical noise model.
% ARG noise : the noise structure for which the log likelihood is required.
% ARG mu : input mean locations for the log likelihood.
% ARG varSigma : input variance locations for the log likelihood.
% ARG y : target locations for the log likelihood.
%
% SEEALSO : orderedNoiseParamInit, orderedNoiseLikelihood, noiseLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


D = size(y, 2);

L = 0;
fact = sqrt(2)/2;
c = 1./sqrt(noise.variance + varsigma);
for j = 1:D
  % Do lowest category first
  index = find(y(:, j) == 0);
  if ~isempty(index)
    mu(index, j) = mu(index, j) + noise.bias(j) ;
    L = L + sum(lnCumGaussian(-mu(index, j).*c(index, j)));
  end
  
  % Intermediate categories
  index = find(y(:, j) > 0 & y(:, j) < noise.C-1);
  if ~isempty(index)
    for i = index'
      mu(i, j) = mu(i, j) + noise.bias(j) - sum(noise.widths(1:y(i, j)-1));
      u = mu(i, j).*c(i, j);
      uprime = (mu(i, j) - noise.widths(y(i, j))).* c(i, j); 
      L = L + lnDiffCumGaussian(u, uprime); 
    end
  end
  % Highest category
  index = find(y(:, j) == noise.C-1);
  if ~isempty(index)
    for i = index'
      mu(i, j) = mu(i, j) + noise.bias(j) - sum(noise.widths(1:y(i, j)-1));
    end
    L = L + sum(lnCumGaussian(mu(index, j).*c(index, j)));
  end
end
  
