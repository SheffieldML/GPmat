function options = fgplvmOptions(varargin);

% FGPLVMOPTIONS Return default options for FGPLVM model.
% FORMAT
% DESC returns the default options in a structure for a FGPLVM model.
% ARG approx : approximation type, either 'ftc' (no approximation),
% 'dtc' (deterministic training conditional), 'dtcvar', variational
% sparse approximation, 'fitc' (fully
% independent training conditional) or 'pitc' (partially
% independent training conditional.
% RETURN options : structure containing the default options for the
% given approximation type.
%
% SEEALSO : fgplvmCreate
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2009

% FGPLVM

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
