function options = gpOptions(approx);

% GPOPTIONS Return default options for GP model.

% FGPLVM

if nargin < 1
  options.approx = 'ftc';
else
  options.approx = approx;
end

% Select type of optimiser.
options.optimiser = 'conjgrad';

% Set to 1 to learn output scales.
options.learnScales = 0;

% Set to 1 if output processes have a shared variance.
options.isSpherical = 1;

% Set to 1 if there is data missing in the target matrix.
options.isMissingData = 0;

switch options.approx
 case 'ftc'
  options.kern = {'rbf', 'bias', 'white'};
  options.numActive = 0;
  options.beta = [];
 case {'fitc', 'pitc', 'dtc'}
  options.kern = {'rbf', 'bias', 'white'};
  options.numActive = 100;
  options.beta = 1e3;
end

