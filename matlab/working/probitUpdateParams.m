function model = probitUpdateParams(model, index)

% PROBITUPDATEPARAMS Update parameters for probit noise model.

% IVM

if nargin < 2
  index = [1:length(model.diagA)]';
end
sqrtOnePlusA = sqrt(1+model.diagA(index));
model.z(index) = model.y(index).*(model.h(index) + model.bias)./sqrtOnePlusA;
model.alpha(index) = model.y(index).*ngaussian(model.z(index))./(cummGaussian(model.z(index)).*sqrtOnePlusA);
model.nu(index) = model.alpha(index).*(model.alpha(index) + (model.h(index) + model.bias)./(1 + model.diagA(index)));

