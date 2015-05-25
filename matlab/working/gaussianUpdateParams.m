function model = gaussianUpdateParams(model, index)

% GAUSSIANUPDATEPARAMS Update parameters for probit noise model.

% IVM

if nargin < 2
  index = 1:length(model.diagA);
end

model.nu(index) = model.sitePrecision(index)./(1+model.sitePrecision(index).*model.diagA(index));