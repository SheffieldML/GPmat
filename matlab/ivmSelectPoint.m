function [indexSelect, infoChange] = ivmSelectPoint(model, add);

% IVMSELECTPOINT Choose a point for inclusion or removal.

% IVM

if nargin < 2
  % If add is 1, then we are including a point.
  add = 1;
end

switch model.selectionCriterion
 case 'random'
  if add
    indexSelect = ceil(rand(1)*length(model.J));
    infoChange = -.5*sum(log2(1-model.varSigma(indexSelect, :)* ...
                              model.nu(indexSelect, :)), 2);
  else
    indexSelect = ceil(rand(1)*length(model.I));
    infoChange = -.5*sum(log2(1-model.varSigma(indexSelect, :)* ...
                          model.beta(indexSelect, :)+1e-300), 2);
  end
 case 'entropy' 
  delta = ivmComputeInfoChange(model, add);
  [infoChange, indexSelect] = max(delta);
  numSelect = sum(delta==infoChange);
  if numSelect>1
    index1 = find(delta==infoChange);
    index1Select = ceil(rand(1)*numSelect);
    indexSelect = index1(index1Select);
  end
 case 'rentropy' 
  % entropy with first point random
  if length(model.I)
    % if point is already selected select another.
    delta = ivmComputeInfoChange(model, add);
    [infoChange, indexSelect] = max(delta);
  else
    % otherwise select one randomly
    if add
      indexSelect = ceil(rand(1)*length(model.J));
      infoChange = -.5*sum(log2(1-model.varSigma(indexSelect, :)* ...
                                model.nu(indexSelect, :)), 2);
    else
      indexSelect = ceil(rand(1)*length(model.I));
      infoChange = -.5*sum(log2(1-model.varSigma(indexSelect, :)* ...
                                model.beta(indexSelect, :)+1e-300), 2);
    end
  end
end
