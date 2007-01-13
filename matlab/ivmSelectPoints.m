function model = ivmSelectPoints(model, display)

% IVMSELECTPOINTS Selects the point for an IVM.
% FORMAT
% DESC is an alternative selection approach to ivmOptimiseIVM. It
% performs EP style updates and point swaps as well as point
% selection.
% ARG model : model for which the points are being selected.
% ARG display : whether or not to display the details of the
% optimisaiton (set to 0, 1 or 2, default value is 1).
% RETURN model : model structure with the active set selected.
%
% SEEALSO : ivmOptimiseIVM, ivmSelectPoint, ivmAddPoint,
% ivmSelectVisualise, ivmRemovePoint, ivmEpUpdatePoint
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% IVM

tol = 1e-3;
if nargin < 2
  display = 1;
end

model = ivmInit(model);
numData = size(model.X, 1);
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
if length(model.J)>1 & ~strcmp(model.selectionCriterion, 'random')
  deltaInfo = tol+1;
  while abs(deltaInfo) > tol  
    [indexSelectRemove, infoChangeRemove] = ivmSelectPoint(model, 0);
    dataIndexSelectRemove = model.I(indexSelectRemove);
    model = ivmRemovePoint(model, dataIndexSelectRemove);
    [indexSelect, infoChangeAdd] = ivmSelectPoint(model, 1);  
    dataIndexSelectAdd = model.J(indexSelect);
    model = ivmAddPoint(model, dataIndexSelectAdd);
    deltaInfo = infoChangeRemove + infoChangeAdd;
    fprintf('Swapping inclusion %d: point %d for point %d, change %2.4f.\n', ...
            indexSelectRemove, ...
            dataIndexSelectRemove, dataIndexSelectAdd, deltaInfo);
  end   
end
lengthIndex = length(model.I);
betaChange = 1;
oldBeta = model.beta;
counter = 0;
while betaChange > tol
  I = model.I;
  counter = counter + 1;
  for i = I'
    model = ivmEpUpdatePoint(model, i);
    if ~rem(counter, 101)
      model = ivmComputeLandM(model);  
      [model.mu, model.varSigma] = ivmPosteriorMeanVar(model, model.X);
      model = ivmUpdateNuG(model);
    end
  end
  model = ivmComputeLandM(model);  
  [model.mu, model.varSigma] = ivmPosteriorMeanVar(model, model.X);
  model = ivmUpdateNuG(model);
  betaChange = full(max(abs(model.beta - oldBeta)));
  fprintf('EP Update %d, beta change %2.4f\n', counter, betaChange)
  oldBeta =  model.beta;
end
