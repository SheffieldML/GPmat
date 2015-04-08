function model = updateCholesky(model, index, activePoint, inactiveSet)

% UPDATECHOLESKY Update Cholesky factor L and matrix M.

% IVM

if length(model.activeIndex) > 0
  % Update model.L, model.M, model.h and model.A with selected point
  lvec = sqrt(model.sitePrecision(index))*model.M(:, index);
  l = sqrt(1 + model.sitePrecision(index)*model.diagK(index) - lvec'*lvec);
  model.Kstore(:, activePoint) = kernel(model.X, model.lntheta, model.kernelType, ...
					model.X(index, :));
  % introduce diagonal term
  model.Kstore(index, activePoint) = model.Kstore(index, activePoint) + exp(model.lntheta(3)); 
  % muvec = 1/l*(sqrt(model.sitePrecision(index))*model.Kstore(:, activePoint) - model.M'*lvec); % Matthias's
  muvec = (1/l)*sqrt(model.sitePrecision(index))*(model.Kstore(:, activePoint) - model.M'*model.M(:, index)); % Neil's
  
  model.L = [model.L zeros(activePoint-1, 1); lvec' l];
  %  model.M = [model.M; zeros(1, numData)];%
  %  model.M(:, inactiveSet) = [model.M(1:end-1, inactiveSet); muvec(inactiveSet)'];
  model.M = [model.M; muvec'];
  model.diagA(inactiveSet) = model.diagA(inactiveSet) - muvec(inactiveSet).*muvec(inactiveSet);
else
  model.L = sqrt(1 + model.sitePrecision(index)*model.diagK(index));
  model.Kstore(:, 1) = kernel(model.X, model.lntheta, model.kernelType, model.X(index, :));
  
  % It is not always clear that we should add this here.
  model.Kstore(index, 1) = model.Kstore(index, 1) + exp(model.lntheta(3));
  
  model.M = (1/model.L*(sqrt(model.sitePrecision(index))*model.Kstore(:, 1)))';
  model.diagA = model.diagA - model.M'.*model.M';
end
