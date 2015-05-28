function [indexSelect, infoChange] = ivmSelectPoint(model, add);

% IVMSELECTPOINT Choose a point for inclusion or removal.
% FORMAT
% DESC identifies the next point for inclusion or removal.
% ARG model : IVM structure for which the next point is being
% selected.
% ARG add : flag which indicates whether or not we are adding a
% point. If we are not adding we are assumed to be removing a point
% (default is true).
%
% SEEALSO : ivmOptimiseIvm, ivmSelectPoints, ivmCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2007

% IVM

if nargin < 2
  % If add is 1, then we are including a point.
  add = true;
end

switch model.selectionCriterion
 case 'random'
  if add
    indexSelect = ceil(rand(1)*length(model.J));
    infoChange = -.5*sum(log2(1-model.varSigma(indexSelect, :).* ...
                              model.nu(indexSelect, :)), 2);
  else
    indexSelect = ceil(rand(1)*length(model.I));
    infoChange = -.5*sum(log2(1-model.varSigma(indexSelect, :).* ...
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
