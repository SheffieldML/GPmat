function options = gpOptions(approx);

% GPOPTIONS Return default options for GP model.
% FORMAT
% DESC returns the default options in a structure for a GP model.
% ARG approx : approximation type, either 'ftc' (no approximation),
% 'dtcvar' (variational sparse approximation) 'dtc' (deterministic training
% conditional), 'fitc' (fully independent training conditional) or 'pitc'
% (partially independent training conditional).
% RETURN options : structure containing the default options for the
% given approximation type.
%
% SEEALSO : gpCreate
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2007, 2009

% GP

if nargin < 1
  options.approx = 'ftc';
else
  options.approx = approx;
end

% Select type of optimiser.
options.optimiser = optimiDefaultOptimiser;

% Set to true to learn output scales.
options.learnScales = false;

% Set to true to scale outputs to variance 1.
options.scale2var1 = false;

% Set to true to optimise beta.
switch approx
 case 'ftc'
  options.optimiseBeta = false;
 otherwise
  options.optimiseBeta = true;
end

% Set to a given mean function to have a mean function.
options.meanFunction = [];
% Options structure for mean function options.
options.meanFunctionOptions = [];

% Set to 1 if output processes have a shared variance.
options.isSpherical = 1;

% Set to 1 if there is data missing in the target matrix.
options.isMissingData = 0;

switch options.approx
 case 'ftc'
  % bog-standard kernel.
  options.kern = {'rbf', 'bias', 'white'};
  options.numActive = 0;
  options.beta = [];
 case {'fitc', 'pitc', 'dtc', 'dtcvar'}
  options.kern = {'rbf', 'bias', 'white'};
  options.numActive = 100;
  options.beta = 1e3;

  % Option to fix the inducing variables to other latent points.
  options.fixInducing = 0;
  options.fixIndices = [];
end

options.computeS = false;
