function options = rbfperiodicOptions(numHidden)

% RBFPERIODICOPTIONS Create a default options structure for the RBFPERIODIC model.
% FORMAT
% DESC creates a default options structure for the periodic radial basis function model
% structure.
% RETURN options : the default options structure.
%
% SEEALSO : rbfperiodicCreate, modelOptions
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

if nargin < 1
  numHidden = 20;
end
options.hiddenDim = numHidden;
options.optimiser = optimiDefaultOptimiser;
