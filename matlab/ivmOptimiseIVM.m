function model = ivmOptimiseIVM(model, display)

% IVMOPTIMISEIVM Selects the points for an IVM model.

% IVM

if nargin < 2
  display = 1;
end
tol = 1e-4;
model = ivmInit(model);
numData = size(model.X, 1);

if model.d > numData
  warning('Active set size is larger than data-set, resetting to data-set size');
  model.d = numData;
end

% Start with all data-points inactive.
model.J = 1:numData;

% Set first infoChange to NaN
infoChange(1) = NaN;
dVal = model.d;
for k = 1:dVal
    
  [indexSelect, infoChange(k)] = ivmSelectPoint(model);
  dataIndexSelect = model.J(indexSelect);
  model = ivmAddPoint(model, dataIndexSelect);
  if display
    ivmSelectVisualise(model, display, k, dataIndexSelect);
  end
end
if model.epUpdate
  lengthIndex = length(model.I);
  betaChange = 1;
  oldBeta = model.beta;
  counter = 0;
  while betaChange > tol
    I = model.I';
    counter = counter + 1;
    for i = I; 
      model = ivmEpUpdatePoint(model, i);
    end
    betaChange = full(max(abs(model.beta - oldBeta)));
    fprintf('EP Update %d, beta change %2.4f\n', counter, betaChange)
    oldBeta =  model.beta;
  end
end