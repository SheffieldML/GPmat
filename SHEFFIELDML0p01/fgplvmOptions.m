function options = fgplvmOptions(varargin);

% FGPLVMOPTIONS Return default options for FGPLVM model.
%
%	Description:
%
%	OPTIONS = FGPLVMOPTIONS(APPROX) returns the default options in a
%	structure for a FGPLVM model.
%	 Returns:
%	  OPTIONS - structure containing the default options for the given
%	   approximation type.
%	 Arguments:
%	  APPROX - approximation type, either 'ftc' (no approximation),
%	   'dtc' (deterministic training conditional), 'dtcvar', variational
%	   sparse approximation, 'fitc' (fully independent training
%	   conditional) or 'pitc' (partially independent training
%	   conditional.
%	
%
%	See also
%	FGPLVMCREATE


%	Copyright (c) 2005, 2009 Neil D. Lawrence


% Get default options from Gaussian process.
options = gpOptions(varargin{:});

% switch optimiser to the OPTIMI specified default.
options.optimiser = optimiDefaultOptimiser;

% How to initialise X.
options.initX = 'ppca';

% What prior on the latent space (only affects starting points if
% dynamics are used.
options.prior = 'gaussian';

% Whether or not to use back constraints
options.back = [];
options.backOptions = [];

% If set to 1 set initial back constraint mapping to match
% initialisation on X. Default is to use standard model initialisation.
options.optimiseInitBack = 0;

% Optional prior on the inducing inputs, X_u
options.inducingPrior = [];
