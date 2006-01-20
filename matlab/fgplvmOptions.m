function options = fgplvmOptions(varargin);

% FGPLVMOPTIONS Return default options for FGPLVM model.

% FGPLVM

options = gpOptions(varargin{:});

options.initX = 'ppca';
options.prior = 'gaussian';
options.back = [];
options.backOptions = [];
options.optimiseInitBack = 1;
options.inducingPrior = 'gaussian';
