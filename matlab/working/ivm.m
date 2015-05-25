function model = ivm(X, y, kernelType, noiseType, selectionCriterion, d)

% IVM Initialise an IVM model.

% IVM

model.X = X;
model.y = y;
model.siteMean = [];
model.sitePrecision = [];
model.nu = zeros(size(y));
model.diagK = zeros(size(y));
model.diagA = zeros(size(y));
model.Kstore = [];
model.activeIndex = [];
model.inactiveIndex = [];
model.L = [];
model.M = [];
model.z = zeros(size(y));
model.h = zeros(size(y));
model.alpha = zeros(size(y));
model.kernelType = kernelType;
model.noiseType = noiseType;
model.selectionCriterion = selectionCriterion;
model.lntheta = initTheta(kernelType, size(X, 2));
model.d = d;

switch noiseType
 
 case 'gaussian'
  model.bias = mean(y);
 
 otherwise
  nClass1 = sum(y==1);
  nClass2 = sum(y==-1);
  model.bias = invCummGaussian(nClass1./(nClass2+nClass1));

end

switch selectionCriterion
 case 'none'
  numData = size(X, 1);
  model.activeIndex = (1:numData);
 otherwise
  
end