function g = probitNoiseGradientParam(noise, mu, varsigma, y)

% PROBITNOISEGRADIENTPARAM Gradient of the probit noise model's parameters.

% IVM

c = y./sqrt(1+ varsigma);
for i = 1:size(mu, 2)
  u(:, i) = c(:, i).*(mu(:, i) + noise.bias(i));
end
denom = cumGaussian(u);
%g = sum(c.*ngaussian(u)./denom, 1);
g = sum(c./(sqrt(2*pi)*(exp(1/2*u.^2)-.5*erfcx(sqrt(2)/2*u))));
