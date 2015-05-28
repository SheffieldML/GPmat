function delta = ivmComputeInfoChange(model, add)

% IVMCOMPUTEINFOCHANGE Compute the information change associated with each point.
% FORMAT
% DESC computes the information change associated with adding or
% deleting a point. This function is used to select which points
% are added or deleted.
% ARG model : the model for which the information change is to be
% comptuted.
% ARG add : flag which states whether the points are to be added or
% subtracted. Set to true if points are to be added, false
% otherwise (default value is true).
% RETURN delta : vector of information change values associated
% with the different points.
%
% SEEALSO : ivmCreate, ivmSelectPoint
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2007

% IVM

if nargin < 2
  add = true;
end

if add
  switch model.selectionCriterion
   case {'entropy', 'rentropy'}
    if model.noise.spherical
      % Noise model leads to constant values for beta.
      delta = -.5*size(model.y, 2).*sum(log2(1-model.varSigma(model.J, 1).* ...
                                             model.nu(model.J, 1)+1e-300), 2);
    else
      delta = -.5*sum(log2(1-model.varSigma(model.J, :).* ...
                           model.nu(model.J, :)+1e-300), 2);
    end
   otherwise
    error(['Selection criterion ' model.selectionCriterion ' not yet implemented'])
  end
else
  
  switch model.selectionCriterion
   case {'entropy', 'rentropy'}
    delta = .5*sum(log2(1-model.varSigma(model.I, :).*model.beta(model.I, ...
						  :)+1e-300), 2);
   otherwise
    error(['Selection criterion ' model.selectionCriterion ' not yet implemented'])
    
  end
end
