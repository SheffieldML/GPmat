function options = gpReversibleDynamicsOptions(approx);

% GPREVERSIBLEDYNAMICSOPTIONS Return default options for GP reversible dynamics model.

% FGPLVM


options = gpOptions(approx);
type = {'cmpnd', {'tensor', 'rbf', 'mlp'}, 'white'};
options.kern = kernCreate(4, type);
options.kern.comp{1} = kernSetIndex(options.kern.comp{1}, 1, [1 2]);
options.kern.comp{1} = kernSetIndex(options.kern.comp{1}, 2, [3 4]);
options.kern.comp{1}.comp{1}.inverseWidth = 0.2;
options.kern.comp{1}.comp{1}.variance = 0.001;
options.kern.comp{1}.comp{2}.variance = 2/pi;
options.kern.comp{1}.comp{2}.weightVariance = 1000;
options.kern.comp{1}.comp{2}.biasVariance = eps;
switch approx
  case 'ftc'
   options.kern.comp{2}.variance = 1e-6;
 case {'fitc', 'pitc', 'dtc', 'dtcvar'}
  options.kern.comp{2}.variance = 0.5e-6;
  options.beta = 2e6;
end
