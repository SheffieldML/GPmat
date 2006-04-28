function options = gpOptions(approx);

% GPOPTIONS Return default options for GP model.
%
% options = gpOptions(approx);
%

% Copyright (c) 2006 Neil D. Lawrence
% gpOptions.m version 1.3



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
  % bog-standard kernel.
  options.kern = {'rbf', 'bias', 'white'};
  options.numActive = 0;
  options.beta = [];
 case {'fitc', 'pitc', 'dtc'}
  options.kern = {'rbf', 'bias', 'white'};
  options.numActive = 100;
  options.beta = 1e3;

  % Option to fix the inducing variables to other latent points.
  options.fixInducing = 0;
  options.fixIndices = [];
    
 case 'nftc'
     % Like FTC but with presition noise beta
  options.kern = {'rbf', 'bias', 'white'};
  options.numActive = 0;
  options.beta = 1e3;     
  
end

options.KLCorrectionTerm = 0; % A setting for the KL correction
