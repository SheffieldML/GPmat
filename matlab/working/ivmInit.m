function model = ivmInit(model, d)

% IVMINIT Initialise the IVM model.

% IVM

if nargin < 2
  d = model.d;
end

% Get number of data.
numData = size(model.y, 1);

% Initialise kernel storage.
model.Kstore = zeros(numData, d);

% Initialise Indices
switch model.selectionCriterion
 case 'none'
  model.activeIndex = 1:size(model.X, 1);
 otherwise
  model.activeIndex = [];
end


% Initialise parameters
model.diagK = kerneldiag(model.X, model.lntheta, model.kernelType);

% Initialise site precision and mean
switch model.noiseType  
 
 case 'probit'
  model.diagA = model.diagK;
  model.siteMean = sparse(zeros(numData, 1));
  model.sitePrecision = sparse(zeros(numData, 1));
  model.h = zeros(numData, 1);
  model.z = zeros(numData, 1);
  model.alpha = zeros(numData, 1);
  model.nu = zeros(numData, 1);
  model = probitUpdateParams(model);
 case 'gaussian'
  model.diagA = model.diagK;
  model.siteMean = model.y;
  model.sitePrecision = ones(size(model.siteMean))*1/exp(model.lntheta(3));
  model = gaussianUpdateParams(model);

end

