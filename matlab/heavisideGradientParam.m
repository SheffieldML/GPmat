function g = heavisideGradientParam(model, params)

% HEAVISIDEGRADIENTPARAM Gradient of the heaviside noise model's parameters.

% IVM

c = model.y./sqrt(model.varSigma);
for i = 1:size(model.mu, 2)
  u(:, i) = c(:, i).*(model.mu(:, i) + model.noise.bias(i));
end
denom = ((1-2*model.noise.eta)*cumGaussian(u)+model.noise.eta);
gnoise.bias = (1-2*model.noise.eta).*sum(c.*ngaussian(u)./denom, 1);
gnoise.eta = sum((-2*cumGaussian(u)+1)./denom, 1);
gnoise.eta = gnoise.eta.*(model.noise.eta.*(1-model.noise.eta));
g = [gnoise.eta gnoise.bias];
