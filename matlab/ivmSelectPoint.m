function   [indexSelect, infoChange] = ivmSelectPoint(model);

% IVMSELECTPOINT Choose a point for inclusion from the inactive set.

% IVM

switch model.selectionCriterion
 case 'random'
  indexSelect = ceil(rand(1)*length(model.J));
  infoChange = -.5*log2(1-model.varSigma(indexSelect)*model.nu(indexSelect));
 case 'entropy' 
  delta = computeInfoChange(model);
  [infoChange, indexSelect] = max(delta);
  if sum(delta==infoChange)==length(delta);
    indexSelect = ceil(rand(1)*length(delta));
  end
end
