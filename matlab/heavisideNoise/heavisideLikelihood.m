function L = heavisideLikelihood(noise, mu, varsigma, y)

% HEAVISIDELIKELIHOOD Likelihood of data under heaviside noise model.

% IVM

D = size(y, 2);
L = zeros(size(mu));
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
  L(:, i) = (1-2*noise.eta)...
            *cumGaussian((y(:, i).*mu(:, i))...
                         ./(sqrt(varsigma(:, i))))+noise.eta;
end
