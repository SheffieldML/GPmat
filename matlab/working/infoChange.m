function delta = infoChange(model)

% INFOCHANGE Compute the information change associated with each point.

switch model.selectionCriterion
 case 'entropy'

  delta = -.5*log2(1-model.diagA(model.inactiveIndex).* ...
		   model.nu(model.inactiveIndex));
 otherwise
  error(['Selection criterion ' model.selectionCriterion [' not yet' ...
		    ' implemented'])
end