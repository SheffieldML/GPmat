function [nu, g] = probitNoiseUpdateParams(noise, mu, varsigma, y, index)

% PROBITNOISEUPDATEPARAMS Update parameters for probit noise model.

% IVM

c = y(index, :)./sqrt(1+varsigma(index, :));

counter = 0;
u = zeros(size(c));
for i = index
  counter = counter + 1;
  u(counter, :) = c(counter, :).*(mu(i, :) + noise.bias);
end

% This version handles large -ve values
% Note, erfcx(x) = exp(x^2) - exp(x^2)erf(x);
g = c./(sqrt(2*pi)*(exp(1/2*u.^2)-.5*erfcx(sqrt(2)/2*u)));

% This version is normal
%g(index, :) = c(index, :).*ngaussian(u(index, :))./(cumGaussian(u(index, :)));
nu = g.*(g + c.*u);

