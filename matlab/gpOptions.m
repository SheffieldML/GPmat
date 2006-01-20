function options = gpOptions(approx);

% GPOPTIONS Return default options for FGPLVM model.

% FGPLVM

if nargin < 1
  options.approx = 'ftc';
else
  options.approx = approx;
end

options.optimiser = 'conjgrad';
options.learnScales = 0;

switch options.approx
 case 'ftc'
  options.kern = {'rbf', 'bias', 'white'};
  options.numActive = 0;
  options.beta = [];
 case {'fitc', 'pitc', 'dtc'}
  options.kern = {'rbf', 'bias', 'white'};
  options.numActive = 100;
  options.beta = 1e-3;
end

