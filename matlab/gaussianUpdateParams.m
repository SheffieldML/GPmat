function model = gaussianUpdateParams(model, index)

% GAUSSIANUPDATEPARAMS Update parameters for probit noise model.

% IVM

D = size(model.y, 2);
model.nu(index, :) = 1./(model.noise.sigma2+model.varSigma(index, :));
for i = 1:D
  model.g(index, i) = model.y(index, i) - model.mu(index, i) - ...
      model.noise.bias(i);
end
model.g(index, :) = model.g(index, :).*model.nu(index, :);