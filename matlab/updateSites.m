function model = updateSites(model, index)

% UPDATESITES Update site parameters.

% IVM

%/~
if any(1./model.nu(index, :)<model.varSigma(index, :))
  warning('nu^-1 is less than varsigma')
end
%~/
model.beta(index, :) = 1./(1./model.nu(index, :)-model.varSigma(index, :));
model.beta(index, :) = thetaConstrain(model.beta(index, :));
model.m(index, :) = model.mu(index, :) + model.g(index, :)./model.nu(index, :);
