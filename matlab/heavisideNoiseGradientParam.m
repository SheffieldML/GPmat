function g = heavisideNoiseGradientParam(noise, mu, varsigma, y)

% HEAVISIDENOISEGRADIENTPARAM Gradient of the heaviside noise model's parameters.

% IVM

c = y./sqrt(varsigma);
denom = zeros(size(c));
u = zeros(size(c));
for i = 1:size(mu, 2)
  u(:, i) = c(:, i).*(mu(:, i) + noise.bias(i));
  denom(:, i) = ((1-2*noise.eta)*cumGaussian(u(:, i))+noise.eta);
end

gnoise.bias = (1-2*noise.eta).*sum(c.*ngaussian(u)./denom, 1);
gnoise.eta = sum(sum((-2*cumGaussian(u)+1)./denom, 1));
gnoise.eta = gnoise.eta.*(noise.eta.*(1-2*noise.eta));
g = [gnoise.eta gnoise.bias];
