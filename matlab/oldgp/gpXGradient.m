function gx = gpXGradient(model)

% GPXGRADIENT Gradient of the kernel with respect to its points.

% GP

xDim = size(model.X, 2);
numData = size(model.X, 1);
%g = zeros(numData, xDim);
gx = kernGradX(model.kern, model.X, model.X);
% factor of 2 accounts for the fact that the covariance is symmetric
gx = gx*2;
% gx has assumed that n is not in model.I, fix that here.
dgx = kernDiagGradX(model.kern, model.X);
for i = 1:numData
  gx(i, :, i) = dgx(i, :);
end
