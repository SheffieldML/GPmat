function options = mlpOptions(numHidden)

% MLPOPTIONS Options for the multi-layered perceptron.
% FORMAT
% DESC returns the default options for a multi-layer perceptron.
% RETURN options : default options structure for Multi-layer
% peceptron.
%
% FORMAT
% returns the default options for a multi-layer perceptron given a
% number of hidden units.
% ARG numHidden : number of hidden units.
% RETURN options : default options structure for Multi-layer
% peceptron with the specified number of hidden units.
%
% SEEALSO : mlpCreate, mlp
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

if nargin < 1
  numHidden = 20;
end
options.hiddenDim = numHidden;
options.activeFunc = 'linear';
options.optimiser = optimiDefaultOptimiser;
