function [nu, g] = heavisideNoiseUpdateParams(noise, mu, varsigma, y, index)

% HEAVISIDENOISEUPDATEPARAMS Update parameters for heaviside noise model.

% IVM

c = y(index, :)./sqrt(varsigma(index, :));
u = zeros(size(c));
for i = 1:size(mu, 2)
  u(:, i) = c(:, i).*(mu(index, i) + noise.bias(i));
  g(:, i) = c(:, i).*ngaussian(u(:, i))...
            ./(cumGaussian(u(:, i))+ ...
               noise.eta(i)/(1-2*noise.eta(i)));
end
nu = g.*(g + c.*u);

