function model = probitUpdateSites(model, index)

% PROBITUPDATESITES Update site parameters for probit model.

% IVM

model = probitUpdateParams(model, index);

model.sitePrecision(index) = model.nu(index)./(1-model.diagA(index).*model.nu(index));
model.siteMean(index) = model.h(index) + model.alpha(index)./model.nu(index);
