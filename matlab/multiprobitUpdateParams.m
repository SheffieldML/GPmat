function model = multiprobitUpdateParams(model, index)

% MULTIPROBITUPDATEPARAMS Update parameters for probit noise model.


model.c(index, :) = model.y(index, :)./sqrt(1+model.varSigma(index, :));
model.u(index, :) = model.c(index, :)...
    .*(model.mu(index, :) + repmat(model.noise.bias, length(index), 1));


% This version handles large -ve values
% Note, erfcx(x) = exp(x^2) - exp(x^2)erf(x);
model.g(index, :) = model.c(index, :)...
    ./(sqrt(2*pi)*(exp(1/2*model.u(index, :).^2)-.5*erfcx(sqrt(2)/2* ...
						  model.u(index, :))));

% This version is normal
%model.g(index, :) = model.c(index, :).*ngaussian(model.u(index, :))./(cummGaussian(model.u(index, :)));
model.nu(index, :) = model.g(index, :).*(model.g(index, :) + model.c(index, :).*model.u(index, :));

