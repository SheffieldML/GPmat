function g = probitGradientParam(model, params)

% PROBITGRADIENTPARAM Gradient of the probit noise model's parameters.

% IVM

c = model.y./sqrt(1+ model.varSigma);
for i = 1:size(model.mu, 2)
  u(:, i) = c(:, i).*(model.mu(:, i) + model.noise.bias(i));
end
denom = cumGaussian(u);
%g = sum(c.*ngaussian(u)./denom, 1);
g = sum(c./(sqrt(2*pi)*(exp(1/2*u.^2)-.5*erfcx(sqrt(2)/2*u))));
