function model = probitUpdateParams(model, index)

% PROBITUPDATEPARAMS Update parameters for probit noise model.

% IVM

model.c(index, :) = model.y(index, :)./sqrt(1+model.varSigma(index, :));

for i = index
  model.u(i, :) = model.c(i, :).*(model.mu(i, :) + model.noise.bias);
end

% This version handles large -ve values
% Note, erfcx(x) = exp(x^2) - exp(x^2)erf(x);
model.g(index, :) = model.c(index, :)./(sqrt(2*pi)*(exp(1/2*model.u(index, :).^2)-.5*erfcx(sqrt(2)/2*model.u(index, :))));

% This version is normal
%model.g(index, :) = model.c(index, :).*ngaussian(model.u(index, :))./(cumGaussian(model.u(index, :)));
model.nu(index, :) = model.g(index, :).*(model.g(index, :) + model.c(index, :).*model.u(index, :));

