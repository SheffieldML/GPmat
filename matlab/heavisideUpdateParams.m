function model = heavisideUpdateParams(model, index)

% HEAVISIDEUPDATEPARAMS Update parameters for heaviside noise model.

% IVM

model.c(index, :) = model.y(index, :)./sqrt(model.varSigma(index, :));
model.u(index, :) = model.c(index, :).*(model.mu(index, :) + model.noise.bias);
model.g(index, :) = model.c(index, :).*ngaussian(model.u(index, :))./(cumGaussian(model.u(index, :))+model.noise.eta/(1-2*model.noise.eta));
model.nu(index, :) = model.g(index, :).*(model.g(index, :) + model.c(index, :).*model.u(index, :));

