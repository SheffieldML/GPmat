function options = rbfperiodicOptions(numHidden)

% RBFPERIODICOPTIONS Create a default options structure for the RBFPERIODIC model.
%
%	Description:
%
%	OPTIONS = RBFPERIODICOPTIONS creates a default options structure for
%	the periodic radial basis function model structure.
%	 Returns:
%	  OPTIONS - the default options structure.
%	
%
%	See also
%	RBFPERIODICCREATE, MODELOPTIONS


%	Copyright (c) 2007 Neil D. Lawrence


if nargin < 1
  numHidden = 20;
end
options.hiddenDim = numHidden;
options.optimiser = optimiDefaultOptimiser;